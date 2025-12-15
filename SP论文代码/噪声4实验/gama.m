clear;close all;clc;
N=5;
L=5000;
Begin=500;%拟收敛的噪声样本点数
Length=500;%EM拟估计的噪声样本点数
UU=randn(N,L);
w0=randn(N,1);
orders=2;%建模为二混合高斯
% MCC的参数
sigma_MCC=10;
% MEE的参数
MEE_Length=10;
sigma_MEE=10;
   %% 各种步长的设置

    mu_Rand=0.21;%随机梯度的步长
    mu_MEE=0.6;%MEE的步长
    mu_LMS=0.005;%普通LMS的步长
    mu_MCC=0.006;%MCC的步长

  tic
  for kk=1:10%做多次实验取平均值    
%%  噪声
    VV= gamrnd(8,4,1,L);% 伽马噪声(alpha=8,beta=4)
    VV=VV-mean(VV);
 %%  噪声构建完毕
    DD=w0'*UU+VV;
    %% EM算法拟估计噪声模型均值，标准差和权重
    [MIU_EM,SIGMA_EM,WEIGHT_EM]=EMtoMixGaussModel(VV(1:Length),orders);
    SIGMA_EM=sqrt(SIGMA_EM);
    R2=0;
for k=1:orders
    R2=R2+WEIGHT_EM(k)/SIGMA_EM(k)^2;
end
    %% 预先定义均值、方差和权重变量
    MIU_Rand=MIU_EM;
    SIGMA_Rand=SIGMA_EM;
    WEIGHT_Rand=WEIGHT_EM;

    %% 利用LMS进行拟收敛
    w_initial=randn(N,1);
    w_Rand=w_initial;
    w_LMS=w_initial;
    w_MCC=w_initial;
    w_MEE=w_initial;
        %% 正式迭代
 for i=1:L
        dk=DD(i);
        uk=UU(:,i);
       %% LMS
        Err_LMS(kk,i)=(w_LMS-w0)'*(w_LMS-w0);
        ek_LMS=dk-w_LMS'*uk;
        w_LMS=w_LMS+mu_LMS*ek_LMS*uk;
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
            w_Rand=w_Rand+mu_Rand*R1*uk;
           end
   %% 对比的算法
   %% MCC
        Err_MCC(kk,i)=(w_MCC-w0)'*(w_MCC-w0);
        ek_MCC = dk- w_MCC'*uk;
          if i<Begin
            w_MCC=w_MCC+mu_LMS*ek_MCC*uk;
        else
        w_MCC = w_MCC + mu_MCC * exp(-(ek_MCC)^2/(2*sigma_MCC^2))*(ek_MCC)*uk;
        end
  %% MEE
        Err_MEE(kk,i)=(w_MEE-w0)'*(w_MEE-w0);
        ek_MEE1 = dk- w_MEE'*uk;
         if i<Begin
            w_MEE=w_MEE+mu_LMS*ek_MEE1*uk;
         else
        for j = 1:MEE_Length
            ek_MEE(j) = DD(i - MEE_Length + j) - w_MEE' * UU(:,i - MEE_Length + j);
        end 
        for ii = 1 : MEE_Length
            for j = 1 : MEE_Length
                uij = UU(:,i - MEE_Length + ii) - UU(:,i - MEE_Length + j);
                eij = ek_MEE(ii) - ek_MEE(j);
                tmp = 1/(MEE_Length^2*sigma_MEE^2) * exp(-(eij^2)/2/sigma_MEE^2) * eij * uij;
                w_MEE = w_MEE + mu_MEE * tmp;
            end           
        end
         end
    TH_Rand(kk,i)=mu_Rand^2*5*R2/(1-(1-mu_Rand*R2)^2);
 end
     %% 输出测验次数
      disp(kk)
  end
  toc
%% 画图
% 算法比较的图示
figure,hold on;
plot(10*log10(mean(Err_LMS)),'-g');
plot(10*log10(mean(Err_MCC)),'-r');  
plot(10*log10(mean(Err_MEE)),'-c');
plot(10*log10(mean(Err_Rand)),'-b');
plot(10*log10(mean(TH_Rand)),'--b');
legend('LMS','MCC (\sigma=10)','MEE (L=10,\sigma=10)','BPD','TH-BPD');
xlabel('Iteration');
ylabel('MSD(dB)');
box on;
figure(2)
[f,xi]=ksdensity(VV);
plot(xi,f,'-r');
grid on;


%% 混合高斯噪声
figure(3)
    Mixedgaussian=zeros(1,L);
    for n=1:L
        Epsilon=rand;
        if Epsilon<WEIGHT_EM(1)
            Mixedgaussian(n)=normrnd(MIU_EM(1),SIGMA_EM(1));
        else
           Mixedgaussian(n)=normrnd(MIU_EM(2),SIGMA_EM(2));
        end
    end
      [f1,xi1]=ksdensity(Mixedgaussian);
      plot(xi1,f1,'-r');
grid on;

