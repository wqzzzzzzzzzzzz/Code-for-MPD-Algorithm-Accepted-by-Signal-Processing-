clear;close all;clc;
N=5;
L=5000;
Begin=500;
UU=randn(N,L);
w0=randn(N,1);
Length=500;%EM拟估计的噪声样本点数
orders=2;
%% 步长的设置
mu1=0.012;%步长
mu2=0.015;
mu3=0.02;
mu4=0.03;

mu_LMS=0.005;%普通LMS的步长
%%  噪声参数的设置
% 均值
miu_real=[3 12];
% 标准差
sigma_real=[1 2];
% 权重
weight_real=[0.5 0.5];

R=0;
for k=1:orders
    R=R+weight_real(k)/sigma_real(k)^2;
end
tic
for kk=1:100%做多次实验取平均值
    %% 混合高斯噪声
    VV=zeros(1,L);
    for n=1:L
        Epsilon=rand;
        if Epsilon<weight_real(1)
            VV(n)=normrnd(miu_real(1),sigma_real(1));
        elseif weight_real(1)<Epsilon&&Epsilon<(weight_real(1)+weight_real(2))
            VV(n)=normrnd(miu_real(2),sigma_real(2));
        end
    end
       VV= VV-mean(VV);
    %%  噪声构建完毕
    DD=w0'*UU+VV;
    %% EM算法拟估计噪声模型均值，标准差和权重
    [MIU_EM,SIGMA_EM,WEIGHT_EM]=EMtoMixGaussModel(VV(1:Length),orders);
    SIGMA_EM=sqrt(SIGMA_EM);
    %% 预先定义随机梯度和小批量梯度的均值、方差和权重变量
    MIU=MIU_EM;
    SIGMA=SIGMA_EM;
    WEIGHT=WEIGHT_EM;

    w_initial=randn(N,1);
    w_LMS=w_initial;
    w1=w_initial;
    w2=w_initial;
    w3=w_initial;
    w4=w_initial;
  
   for i=1:L
        dk=DD(i);
        uk=UU(:,i);
       %% LMS
        Err_LMS(kk,i)=(w_LMS-w0)'*(w_LMS-w0);
        ek_LMS=dk-w_LMS'*uk;
        w_LMS=w_LMS+mu_LMS*ek_LMS*uk;
       %% 新算法
        Err1(kk,i)=(w1-w0)'*(w1-w0);
        Err2(kk,i)=(w2-w0)'*(w2-w0);
        Err3(kk,i)=(w3-w0)'*(w3-w0);
        Err4(kk,i)=(w4-w0)'*(w4-w0);
    
        ek1=dk-w1'*uk;
        ek2=dk-w2'*uk;
        ek3=dk-w3'*uk;
        ek4=dk-w4'*uk;
 
        if i<Begin
            w1=w1+mu_LMS*ek1*uk;
             w2=w2+mu_LMS*ek2*uk;
              w3=w3+mu_LMS*ek3*uk;
               w4=w4+mu_LMS*ek4*uk;

        else
%%
            for k=1:orders
                P1(i,k)=exp(-1*(ek1-MIU(k))^2/(2*SIGMA(k)^2))/(sqrt(2*pi)*SIGMA(k));
            end
            for k=1:orders
                V1(i,k)=WEIGHT(k)*P1(i,k)/(WEIGHT*P1(i,:)');
            end
            R1=0;
            for k=1:orders
                R1=R1+V1(i,k)*((ek1-MIU(k))/SIGMA(k)^2);
            end
            w1=w1+mu1*R1*uk;
            %%

            for k=1:orders
                P2(i,k)=exp(-1*(ek2-MIU(k))^2/(2*SIGMA(k)^2))/(sqrt(2*pi)*SIGMA(k));
            end
            for k=1:orders
                V2(i,k)=WEIGHT(k)*P1(i,k)/(WEIGHT*P1(i,:)');
            end
            R2=0;
            for k=1:orders
                R2=R2+V2(i,k)*((ek2-MIU(k))/SIGMA(k)^2);
            end
            w2=w2+mu2*R2*uk;
            %%

            for k=1:orders
                P3(i,k)=exp(-1*(ek3-MIU(k))^2/(2*SIGMA(k)^2))/(sqrt(2*pi)*SIGMA(k));
            end
            for k=1:orders
                V3(i,k)=WEIGHT(k)*P3(i,k)/(WEIGHT*P3(i,:)');
            end
            R3=0;
            for k=1:orders
                R3=R3+V3(i,k)*((ek3-MIU(k))/SIGMA(k)^2);
            end
            w3=w3+mu3*R3*uk;
            %%

            for k=1:orders
                P4(i,k)=exp(-1*(ek4-MIU(k))^2/(2*SIGMA(k)^2))/(sqrt(2*pi)*SIGMA(k));
            end
            for k=1:orders
                V4(i,k)=WEIGHT(k)*P4(i,k)/(WEIGHT*P4(i,:)');
            end
            R4=0;
            for k=1:orders
                R4=R4+V4(i,k)*((ek4-MIU(k))/SIGMA(k)^2);
            end
            w4=w4+mu4*R4*uk;
          
        end
          %%
        TH1(kk,i)=mu1^2*5*R/(1-(1-mu1*R)^2);
        TH2(kk,i)=mu2^2*5*R/(1-(1-mu2*R)^2);
        TH3(kk,i)=mu3^2*5*R/(1-(1-mu3*R)^2);
        TH4(kk,i)=mu4^2*5*R/(1-(1-mu4*R)^2);

    end
    disp(kk);
end
toc

figure,hold on;
plot(10*log10(mean(Err_LMS)),'-g');
plot(10*log10(mean(Err1)),'-b');
plot(10*log10(mean(Err2)),'-c');
plot(10*log10(mean(Err3)),'-r');
plot(10*log10(mean(Err4)),'-k');

plot(10*log10(mean(TH1)),'--b');
plot(10*log10(mean(TH2)),'--c');
plot(10*log10(mean(TH3)),'--r');
plot(10*log10(mean(TH4)),'--k');

legend('LMS(\eta=0.005)','BPD1(\eta=0.012)','BPD2(\eta=0.015)','BPD3(\eta=0.02)','BPD4(\eta=0.03)','TH-BPD1','TH-BPD2','TH-BPD3','TH-BPD4');
xlabel('Iteration');
ylabel('MSD(dB)');
box on;




figure(2)
[f,xi]=ksdensity(VV);
plot(xi,f,'-r');
grid on;
xlabel('Magnitude');
ylabel('Possibility');
figure(2)
[f,xi]=ksdensity(VV);
plot(xi,f,'-r');
grid on;
xlabel('Magnitude');
ylabel('Possibility');
title('双峰脉冲噪声');
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
         Mixedgaussian= Mixedgaussian-mean(Mixedgaussian);
      [f1,xi1]=ksdensity(Mixedgaussian);
      plot(xi1,f1,'-r');
grid on;
xlabel('Magnitude');
ylabel('Possibility');
title('EM估计噪声');
