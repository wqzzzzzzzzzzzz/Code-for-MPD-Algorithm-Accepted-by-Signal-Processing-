clear;close all;clc;
N=5;
L=5000;
Begin=500;
UU=randn(N,L);
w0=randn(N,1);
Length=500;%EM拟估计的噪声样本点数
orders=2;
orders1=4;
% MCC的参数
sigma_MCC=2;
% MEE的参数
MEE_Length=10;
sigma_MEE=10;
%% 各种步长的设置

mu_Rand=0.00025;%步长
mu_Rand1=0.00014;%步长
mu_LMS=0.005;%普通LMS的步长
  mu_LMF=0.02;
    mu_MEE=0.3;%MEE的步长
    mu_MCC=0.006;%MCC的步长
%%  噪声参数的设置

tic
for kk=1:10%做多次实验取平均值
    %% 混合高斯噪声
    VV=rand(1,L)-0.5;
         %%  噪声构建完毕
    DD=w0'*UU+VV;
       %% EM1
    [MIU_EM1,SIGMA_EM1,WEIGHT_EM1]=EMtoMixGaussModel_two(VV(1:Length),orders);
    SIGMA_EM1=sqrt(SIGMA_EM1);
    %% EM2
        [MIU_EM2,SIGMA_EM2,WEIGHT_EM2]=EMtoMixGaussModel_four(VV(1:Length),orders1);
    SIGMA_EM2=sqrt(SIGMA_EM2);
    %% R_two
    R_two=0;
for k=1:orders
    R_two=R_two+WEIGHT_EM1(k)/SIGMA_EM1(k)^2;
end
    %% R_four
    R_four=0;
for k=1:orders1
    R_four=R_four+WEIGHT_EM2(k)/SIGMA_EM2(k)^2;
end
    %% 
    MIU_Rand=MIU_EM1;
    SIGMA_Rand=SIGMA_EM1;
    WEIGHT_Rand=WEIGHT_EM1;

       MIU_Rand1=MIU_EM2;
    SIGMA_Rand1=SIGMA_EM2;
    WEIGHT_Rand1=WEIGHT_EM2;

    w_initial=randn(N,1);
    w_LMS=w_initial;
    w_Rand=w_initial;
     w_Rand1=w_initial;
  
      w_LMF=w_initial;
    w_MCC=w_initial;
    w_MEE=w_initial;
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
        %%
        Err_Rand1(kk,i)=(w_Rand1-w0)'*(w_Rand1-w0);
        ek_Rand1=dk-w_Rand1'*uk;
        if i<Begin
            w_Rand1=w_Rand1+mu_LMS*ek_Rand1*uk;
        else

            for k=1:orders1
                P_Rand1(i,k)=exp(-1*(ek_Rand1-MIU_Rand1(k))^2/(2*SIGMA_Rand1(k)^2))/(sqrt(2*pi)*SIGMA_Rand1(k));
            end
            for k=1:orders1
                V_Rand1(i,k)=WEIGHT_Rand1(k)*P_Rand1(i,k)/(WEIGHT_Rand1*P_Rand1(i,:)');
            end
            R1=0;
            for k=1:orders1
                R1=R1+V_Rand1(i,k)*((ek_Rand1-MIU_Rand1(k))/SIGMA_Rand1(k)^2);
            end
            w_Rand1=w_Rand1+mu_Rand1*R1*uk;
           end
   %% 对比的算法
       %% LMF
        Err_LMF(kk,i)=(w_LMF-w0)'*(w_LMF-w0);
        ek_LMF = dk- w_LMF'*uk;
          if i<Begin
            w_LMF=w_LMF+mu_LMS*ek_LMF*uk;
        else
        w_LMF = w_LMF + mu_LMF * ek_LMF^3*uk;
        end
   %% MCC
        Err_MCC(kk,i)=(w_MCC-w0)'*(w_MCC-w0);
        ek_MCC = dk- w_MCC'*uk;
          if i<Begin
            w_MCC=w_MCC+mu_LMS*ek_MCC*uk;
        else
        w_MCC = w_MCC + mu_MCC * exp(-(ek_MCC-mean(VV))^2/(2*sigma_MCC^2))*(ek_MCC-mean(VV))*uk;
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
                 TH_Rand(kk,i)=mu_Rand^2*5*R_two/(1-(1-mu_Rand*R_two)^2);
                  TH_Rand1(kk,i)=mu_Rand1^2*5*R_four/(1-(1-mu_Rand1*R_four)^2);
    end
    disp(kk);
end
toc

figure,hold on;
plot(10*log10(mean(Err_LMS)),'-g');
plot(10*log10(mean(Err_LMF)),'-k'); 
plot(10*log10(mean(Err_MCC)),'-r');  
plot(10*log10(mean(Err_MEE)),'-c');
plot(10*log10(mean(Err_Rand)),'-b');
plot(10*log10(mean(Err_Rand1)),'-m');
plot(10*log10(mean(TH_Rand)),'--b');
plot(10*log10(mean(TH_Rand1)),'--m');
legend('LMS','LMF','MCC (\sigma=2)','MEE (L=10,\sigma=10)','BPD','QPD','TH-BPD','TH-QPD');
xlabel('Iteration');
ylabel('MSD(dB)');
ylim([-36 15]);
xlim([0,L]);
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
        if Epsilon<WEIGHT_EM1(1)
            Mixedgaussian(n)=normrnd(MIU_EM1(1),SIGMA_EM1(1));
        else
           Mixedgaussian(n)=normrnd(MIU_EM1(2),SIGMA_EM1(2));
        end
    end
         Mixedgaussian= Mixedgaussian-mean(Mixedgaussian);
      [f1,xi1]=ksdensity(Mixedgaussian);
      plot(xi1,f1,'-r');
grid on;
xlabel('Magnitude');
ylabel('Possibility');

