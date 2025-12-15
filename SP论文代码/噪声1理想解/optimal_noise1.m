clear;close all;clc;
N=5;
L=15000;
Begin=500;
UU=randn(N,L);
w0=randn(N,1);
Length=500;%EM拟估计的噪声样本点数
orders=2;
%% 各种步长的设置
% 
mu_Rand=0.002;%随机梯度的步长
mu_LMS=0.005;%普通LMS的步长
 
%%  噪声参数的设置
% 均值
miu_real=[-4.5 4.5];
% 标准差
sigma_real=[1 2];
% 权重
weight_real=[0.5 0.5];
R2=0;
for k=1:orders
    R2=R2+weight_real(k)/sigma_real(k)^2;
end

tic
for kk=1:50%做多次实验取平均值
 
    r2=0;
%% 混合高斯噪声
    VV=zeros(1,L);
      p_rand = rand(1,L);
    for n=1:L
        if p_rand(n)<weight_real(1)
            VV(n)=normrnd(miu_real(1),sigma_real(1));
        else
            VV(n)=normrnd(miu_real(2),sigma_real(2));
        end
    end
    VV1=VV;
    %%  噪声构建完毕
    DD=w0'*UU+VV;
    %% 预先定义随机梯度和小批量梯度的均值、方差和权重变量
    MIU_Rand=miu_real ;
    SIGMA_Rand=sigma_real;
    WEIGHT_Rand=weight_real;

    w_initial=randn(N,1);
    w_Rand=w_initial;
    w_wiener = inv(UU*UU')  * UU * DD';

    i1 = 0; i2 =0;
    a1=weight_real(1)/sigma_real(1)^2;
    a2=weight_real(2)/sigma_real(2)^2;

for ii = 1 : L
    if p_rand(ii) > weight_real(1)
        i1 = i1 + 1;
        U1(:, i1) = UU(:, ii);
        D1(i1) = w0' * UU(:, ii) + VV1(ii)-miu_real(2);        
    else
        i2 = i2 + 1;
        U2(:, i2) = UU(:, ii);
        D2(i2) = w0' * UU(:, ii) + VV1(ii)-miu_real(1);   
    end
end

    w_em_th = inv(a2 * U1 * U1' /i1 + a1 * U2 * U2' /i2) * (a2 * U1 * D1' /i1 + a1 * U2 * D2' /i2);
   for i=1:L
        dk=DD(i);
        uk=UU(:,i);
            Err_wiener(kk,i) = (w0 - w_wiener)' * (w0- w_wiener);
                Err_em_th(kk,i) = (w0 - w_em_th)' * (w0 - w_em_th);

       %% 新算法
        Err_Rand(kk,i)=(w_Rand-w0)'*(w_Rand-w0);
        ek_Rand=dk-w_Rand'*uk;
        if i<Begin
                     
            w_Rand=w_Rand+mu_LMS*ek_Rand*uk;
        else
            for k=1:orders
                P_Rand(i,k)=exp(-1*(ek_Rand-MIU_Rand(k))^2/(2*SIGMA_Rand(k)^2))/(sqrt(2*pi)*SIGMA_Rand(k));
            end
            for k=1:orders
                V_Rand(i,k)=WEIGHT_Rand(k)*P_Rand(i,k)/(WEIGHT_Rand*P_Rand(i,:)');
            end
            R1=0;
            for k=1:orders
                R1=R1+V_Rand(i,k)*((ek_Rand-MIU_Rand(k))/SIGMA_Rand(k)^2);
            end
             if i < 2500
            w_Rand=w_Rand+mu_Rand*R1*uk;
            else
                w_Rand=w_Rand+mu_Rand*R1*uk/i*900;
            end
        end
    end
    disp(kk);
end
toc

figure,hold on;
plot(10*log10(mean(Err_Rand)),'-b');
plot(10*log10(mean(Err_wiener)),'-g');
plot(10*log10(mean(Err_em_th)),'-r');
legend('BPD', 'Wiener Solution','Optimal Solution (44)');
xlabel('Iteration');
ylabel('MSD(dB)');
box on;
figure(2)
[f,xi]=ksdensity(VV);
plot(xi,f,'-r');
grid on;
xlabel('Magnitude');
ylabel('Possibility');

%% 混合高斯噪声
figure(3)
    Mixedgaussian=zeros(1,L);
    for n=1:L
        Epsilon=rand;
        if Epsilon<weight_real(1)
            Mixedgaussian(n)=normrnd(miu_real(1),sigma_real(1));
        else
           Mixedgaussian(n)=normrnd(miu_real(2),sigma_real(2));
        end
    end
         Mixedgaussian= Mixedgaussian-mean(Mixedgaussian);
      [f1,xi1]=ksdensity(Mixedgaussian);
      plot(xi1,f1,'-r');
grid on;
xlabel('Magnitude');
ylabel('Possibility');



