%%  
L=500;
orders=2;
orders1=4;
%%
subplot(321)
miu_real3=[3 12];
sigma_real3=[1 2];
weight_real3=[0.5 0.5];

%% 混合高斯噪声
    VV=zeros(1,L);
    for n=1:L
        Epsilon=rand;
        if Epsilon<weight_real3(1)
            VV(n)=normrnd(miu_real3(1),sigma_real3(1));
        else
            VV(n)=normrnd(miu_real3(2),sigma_real3(2));
        end
    end
         %% EM3
    [MIU_EM3,SIGMA_EM3,WEIGHT_EM3]=EMtoMixGaussModel3(VV,orders);
    SIGMA_EM3=sqrt(SIGMA_EM3);
% [f,xi]=ksdensity(VV);
xi=-15:0.01:20;
f=(xi<20&xi>-15)*0.5.*(exp(-(xi-3).^2/2)/(sqrt(2*pi))+exp(-(xi-12).^2/2/4)/(sqrt(2*pi))/2);
plot(xi,f,'-r');hold on;
grid on;
%% 
    Mixedgaussian=zeros(1,L);
    for n=1:L
        Epsilon=rand;
        if Epsilon<WEIGHT_EM3(1)
            Mixedgaussian(n)=normrnd(MIU_EM3(1),SIGMA_EM3(1));
        else
           Mixedgaussian(n)=normrnd(MIU_EM3(2),SIGMA_EM3(2));
        end
    end
      [f1,xi1]=ksdensity(Mixedgaussian);
      plot(xi1,f1,'-g');
grid on;
xlim([-5,20]);
xlabel('(1)');
%%
subplot(322)
miu_real3=[3 6];
sigma_real3=[1 2];
weight_real3=[0.5 0.5];

%% 混合高斯噪声
    VV=zeros(1,L);
    for n=1:L
        Epsilon=rand;
        if Epsilon<weight_real3(1)
            VV(n)=normrnd(miu_real3(1),sigma_real3(1));
        else
            VV(n)=normrnd(miu_real3(2),sigma_real3(2));
        end
    end
         %% EM3
    [MIU_EM3,SIGMA_EM3,WEIGHT_EM3]=EMtoMixGaussModel3(VV,orders);
    SIGMA_EM3=sqrt(SIGMA_EM3);
% [f,xi]=ksdensity(VV);
xi=-10:0.01:15;
f=(xi<15&xi>-10)*0.5.*(exp(-(xi-3).^2/2)/(sqrt(2*pi))+exp(-(xi-6).^2/2/4)/(sqrt(2*pi))/2);
plot(xi,f,'-r');hold on;
grid on;
%% 
    Mixedgaussian=zeros(1,L);
    for n=1:L
        Epsilon=rand;
        if Epsilon<WEIGHT_EM3(1)
            Mixedgaussian(n)=normrnd(MIU_EM3(1),SIGMA_EM3(1));
        else
           Mixedgaussian(n)=normrnd(MIU_EM3(2),SIGMA_EM3(2));
        end
    end
      [f1,xi1]=ksdensity(Mixedgaussian);
      plot(xi1,f1,'-g');
grid on;
xlim([-5,15]);
xlabel('(2)');

%%
subplot(323)
    %% 瑞利噪声
    VV=raylrnd(8,1,L);
    %% EM5
    [MIU_EM5,SIGMA_EM5,WEIGHT_EM5]=EMtoMixGaussModel5(VV,orders);SIGMA_EM5=sqrt(SIGMA_EM5);
% [f,xi]=ksdensity(VV);
xi=0:0.1:40;
f=(xi<40&xi>0).*xi/64.*exp(-xi.^2/128);
plot(xi,f,'-r');hold on;grid on;xlim([-15,30]);
%% 二混合拟合噪声
    Mixedgaussian=zeros(1,L);
    for n=1:L
        Epsilon=rand;
        if Epsilon<WEIGHT_EM5(1)
            Mixedgaussian(n)=normrnd(MIU_EM5(1),SIGMA_EM5(1));
        else
           Mixedgaussian(n)=normrnd(MIU_EM5(2),SIGMA_EM5(2));
        end
    end
      [f1,xi1]=ksdensity(Mixedgaussian);plot(xi1,f1,'-g');
grid on;xlim([-10,40]);xlabel('(3)');
%%
subplot(324)
%%  噪声
    VV= gamrnd(8,4,1,L);% 伽玛噪声
    %% EM算法拟估计噪声模型均值，标准差和权重
    [MIU_EM6,SIGMA_EM6,WEIGHT_EM6]=EMtoMixGaussModel6(VV,orders);
    SIGMA_EM6=sqrt(SIGMA_EM6);
% [f,xi]=ksdensity(VV);
xi=0:0.1:70;
f=(xi<70&xi>0).*(1/5040/4^8).*xi.^7.*exp(-xi./4);
plot(xi,f,'-r');
grid on;hold on;
%% Estimate 噪声
    Mixedgaussian=zeros(1,L);
    for n=1:L
        Epsilon=rand;
        if Epsilon<WEIGHT_EM6(1)
            Mixedgaussian(n)=normrnd(MIU_EM6(1),SIGMA_EM6(1));
        else
           Mixedgaussian(n)=normrnd(MIU_EM6(2),SIGMA_EM6(2));
        end
    end
      [f1,xi1]=ksdensity(Mixedgaussian);
    
      plot(xi1,f1,'-g');
grid on;
xlim([-10,85]);
xlabel('(4)');
%%
 subplot(325)
    %% 均匀分布噪声
 VV=rand(1,L)-0.5;
       %% EM71
    [MIU_EM71,SIGMA_EM71,WEIGHT_EM71]=EMtoMixGaussModel71(VV,orders);
    SIGMA_EM71=sqrt(SIGMA_EM71);
    %% EM72
  [MIU_EM72,SIGMA_EM72,WEIGHT_EM72]=EMtoMixGaussModel72(VV,orders1);
    SIGMA_EM72=sqrt(SIGMA_EM72);
    %% Noise
% [f,xi]=ksdensity(VV);
xi=-0.5:0.01:0.5;
f=1*(xi>-0.5&xi<0.5);
plot(xi,f,'-r');  hold on;
grid on;
%% 二混合拟合噪声
    Mixedgaussian=zeros(1,L);
    for n=1:L
        Epsilon=rand;
        if Epsilon<WEIGHT_EM71(1)
            Mixedgaussian(n)=normrnd(MIU_EM71(1),SIGMA_EM71(1));
        else
           Mixedgaussian(n)=normrnd(MIU_EM71(2),SIGMA_EM71(2));
        end
    end
      [f1,xi1]=ksdensity(Mixedgaussian);
      plot(xi1,f1,'-g');grid on;
%% 四混合拟合噪声
 Mixedgaussian=zeros(1,L);
    for n=1:L
        Epsilon=rand;
        if Epsilon<WEIGHT_EM72(1)
            Mixedgaussian(n)=normrnd(MIU_EM72(1),SIGMA_EM72(1));
        elseif WEIGHT_EM72(1)<Epsilon&&Epsilon<(WEIGHT_EM72(1)+WEIGHT_EM72(2))
           Mixedgaussian(n)=normrnd(MIU_EM72(2),SIGMA_EM72(2));
       elseif (WEIGHT_EM72(1)+WEIGHT_EM72(2))<Epsilon&&Epsilon<(WEIGHT_EM72(1)+WEIGHT_EM72(2)+WEIGHT_EM72(3))
           Mixedgaussian(n)=normrnd(MIU_EM72(3),SIGMA_EM72(3));
        else
           Mixedgaussian(n)=normrnd(MIU_EM72(4),SIGMA_EM72(4));
        end
    end
      [f1,xi1]=ksdensity(Mixedgaussian);
      plot(xi1,f1,'-b');
grid on;
xlim([-1,1]);
xlabel('(5)');
%%
subplot(326)
   VV=sin(randn(1,L));
    %% EM81
    [MIU_EM81,SIGMA_EM81,WEIGHT_EM81]=EMtoMixGaussModel81(VV,orders);
    SIGMA_EM81=sqrt(SIGMA_EM81);
    %% EM82
    [MIU_EM82,SIGMA_EM82,WEIGHT_EM82]=EMtoMixGaussModel82(VV,orders1);
    SIGMA_EM82=sqrt(SIGMA_EM82);
%% Noise
L2=100000;
  VV1=sin(randn(1,L2));
[f,xi]=ksdensity(VV1);
% xi3=sin(exp(-(-10:0.01:10).^2/2)/(sqrt(2*pi)));
% f3=xi3.*exp(-xi3.^2/2)/(sqrt(2*pi));
hold on;
plot(xi,f,'-r');
grid on;
%% 拟合二混合噪声
    Mixedgaussian=zeros(1,L);
    for n=1:L
        Epsilon=rand;
        if Epsilon<WEIGHT_EM81(1)
            Mixedgaussian(n)=normrnd(MIU_EM81(1),SIGMA_EM81(1));
        else
           Mixedgaussian(n)=normrnd(MIU_EM81(2),SIGMA_EM81(2));
        end
    end
      [f1,xi1]=ksdensity(Mixedgaussian);
      plot(xi1,f1,'-g');
grid on;
%% 四混合拟合噪声
 Mixedgaussian=zeros(1,L);
    for n=1:L
        Epsilon=rand;
        if Epsilon<WEIGHT_EM82(1)
            Mixedgaussian(n)=normrnd(MIU_EM82(1),SIGMA_EM82(1));
        elseif WEIGHT_EM82(1)<Epsilon&&Epsilon<(WEIGHT_EM82(1)+WEIGHT_EM82(2))
           Mixedgaussian(n)=normrnd(MIU_EM82(2),SIGMA_EM82(2));
       elseif (WEIGHT_EM82(1)+WEIGHT_EM82(2))<Epsilon&&Epsilon<(WEIGHT_EM82(1)+WEIGHT_EM82(2)+WEIGHT_EM82(3))
           Mixedgaussian(n)=normrnd(MIU_EM82(3),SIGMA_EM82(3));
        else
           Mixedgaussian(n)=normrnd(MIU_EM82(4),SIGMA_EM82(4));
        end
    end
      [f1,xi1]=ksdensity(Mixedgaussian);
      plot(xi1,f1,'-b');
grid on;box on;
xlim([-2,2]);
xlabel('(6)');
legend('PDF of Noise','PDF of Estimated Noise by BPD','PDF of Estimated Noise by QPD');
