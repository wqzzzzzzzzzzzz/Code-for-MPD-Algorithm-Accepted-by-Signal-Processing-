clear;close all;clc;
N=5;
L=15000;
Begin=500;%拟收敛的噪声样本点数
Length=500;%EM拟估计的噪声样本点数
UU=randn(N,L);
w0=randn(N,1);
orders=2;%建模为二混合高斯
%% 各种步长的设置
mu_Rand=0.02;%随机梯度的步长
mu_LMS=0.001;%普通LMS的步长
tic
VV= raylrnd(8,1,Begin);% 瑞利噪声(sigma=8)
%% EM算法拟估计噪声模型均值，标准差和权重
[MIU_EM,SIGMA_EM,WEIGHT_EM]=EMtoMixGaussModel(VV(1:Length),orders);
SIGMA_EM=sqrt(SIGMA_EM);
for kk=1:100 %做多次实验取平均值
    %% 噪声
    VV= raylrnd(8,1,L);% 瑞利噪声(sigma=8)
    %% 噪声构建完毕
    DD=w0'*UU+VV;
  
    %% 利用LMS进行拟收敛
    w_initial=randn(N,1);
    w_Rand=w_initial;
    w_wiener = inv(UU*UU') * UU * DD';
    %% 正式迭代
    for i=1:L
        dk=DD(i);
        uk=UU(:,i);
        Err_wiener(kk,i) = (w0 - w_wiener)' * (w0- w_wiener);
        %% 新算法
        Err_Rand(kk,i)=(w_Rand-w0)'*(w_Rand-w0);
        ek_Rand=dk-w_Rand'*uk;
        if i<Begin
            w_Rand=w_Rand+mu_LMS*ek_Rand*uk;
        else
            for k=1:orders
                P_Rand(i,k)=exp(-1*(ek_Rand-MIU_EM(k))^2/(2*SIGMA_EM(k)^2))/(sqrt(2*pi)*SIGMA_EM(k));
            end
            for k=1:orders
                V_Rand(i,k)=WEIGHT_EM(k)*P_Rand(i,k)/(WEIGHT_EM*P_Rand(i,:)');
            end
            R1=0;
            for k=1:orders
                R1=R1+V_Rand(i,k)*((ek_Rand-MIU_EM(k))/SIGMA_EM(k)^2);
            end
            if i < 5000
            w_Rand=w_Rand+mu_Rand*R1*uk;
            else
                w_Rand=w_Rand+mu_Rand*R1*uk/i*2400;
            end
        end
   
    end
    %% 输出测验次数
    %     disp(kk)
    
    U1 = zeros(N, L-5000);
    U2 = U1;
    D1 = zeros(1, L-5000);
    D2 = D1;
    i1 = 0;
    i2 = 0;
    for ii = 5001 : L
        if V_Rand(ii,1) < V_Rand(ii,2)
            i1 = i1 + 1;
            U1(:, i1) = UU(:, ii);
            D1(i1) = w0' * UU(:, ii) + VV(ii)-MIU_EM(2);
        else
            i2 = i2 + 1;
            %%
            U2(:, i2) = UU(:, ii);
            D2(i2) = w0' * UU(:, ii) + VV(ii)-MIU_EM(1);
        end
    end
    a1=WEIGHT_EM(1)/SIGMA_EM(1)^2;
    a2=WEIGHT_EM(2)/SIGMA_EM(2)^2;
    w_em_th = inv(a2 * U1 * U1' /i1 + a1 * U2 * U2' /i2) * (a2 * U1 * D1' /i1 + a1 * U2 * D2' /i2);
    Err_em_th(kk) = (w0 - w_em_th)' * (w0 - w_em_th);
    
end
toc
%% 画图
% 三种算法比较的图示
figure,hold on;
plot(10*log10(mean(Err_Rand)),'-b');
plot(10*log10(mean(Err_wiener)),'-g');
plot(10*log10(mean(Err_em_th))*ones(1,L),'-r');
legend('BPD', 'Wiener Solution','Optimal Solution (44)');
xlabel('Iteration');
ylabel('MSD(dB)');
box on;
