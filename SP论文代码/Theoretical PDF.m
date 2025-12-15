 clear all,close all,clc
%%  
L1=1000;
L2=10000;
L3=100000;
L4=1000000;

figure(1)

%% 均匀分布噪声
 VV1=rand(1,L1)-0.5;
  VV2=rand(1,L2)-0.5;
   VV3=rand(1,L3)-0.5;
    VV4=rand(1,L4)-0.5;

    %% Noise

[f1,xi1]=ksdensity(VV1);
[f2,xi2]=ksdensity(VV2);
[f3,xi3]=ksdensity(VV3);
[f4,xi4]=ksdensity(VV4);
plot(xi1,f1,'-r','linewidth',1);  hold on;
plot(xi2,f2,'-g','linewidth',1);  
plot(xi3,f3,'-k','linewidth',1); 
plot(xi4,f4,'-b','linewidth',1);  
grid on;
xlim([-1,1]);
xlabel('Magnitude');
ylabel('PDF');
legend('1000 points','10000 points','100000 points','1000000 points');
% legend('PDF of Noise 5');

figure(2)

VV1=sin(randn(1,L1));
VV2=sin(randn(1,L2));
VV3=sin(randn(1,L3));
VV4=sin(randn(1,L4));

%% Noise
[f1,xi1]=ksdensity(VV1);
[f2,xi2]=ksdensity(VV2);
[f3,xi3]=ksdensity(VV3);
[f4,xi4]=ksdensity(VV4);

plot(xi1,f1,'-r','linewidth',1);  hold on;
plot(xi2,f2,'-g','linewidth',1);  
plot(xi3,f3,'-k','linewidth',1); 
plot(xi4,f4,'-b','linewidth',1);  
grid on;box on;
xlim([-2,2]);
xlabel('Magnitude');
ylabel('PDF');
legend('1000 points','10000 points','100000 points','1000000 points');
% legend('PDF of Noise 6');
% legend('PDF of Noise','PDF of Estimated Noise by BPD','PDF of Estimated Noise by QPD');