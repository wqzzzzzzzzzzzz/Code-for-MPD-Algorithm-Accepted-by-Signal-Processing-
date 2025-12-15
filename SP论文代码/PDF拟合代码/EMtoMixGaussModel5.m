function [MIU SIGMA WEIGHT]=EMtoMixGaussModel(MixedGaussSeq ,orders)
%% 不同混合高斯噪声所需的初始值设定
% 二混合高斯噪声
orders=2;
MIU=[7 20];%一开始随机给的均值
SIGMA=[9 36];%一开始随机给的方差，注意是方差
WEIGHT=[0.5 0.5];%一开始随机给的权重
% 三混合高斯噪声
% orders=3;
% MIU=[-6 -2 10];%一开始随机给的均值
% SIGMA=[1 3 5];%一开始随机给的方差
% WEIGHT=[0.4 0.4 0.2];%一开始随机给的权重
%% 正式进行EM估计
L=length(MixedGaussSeq);%信号长度
delta=1;%delta的初始条件
P=zeros(L,orders);
while delta>0.0000000001 %结束迭代的终止条件
    %E step
    for i=1:L
        for k=1:orders
            P(i,k)=exp(-1*(MixedGaussSeq(i)-MIU(k))^2/(2*SIGMA(k)))/sqrt(2*pi*SIGMA(k));%概率密度函数
        end
    end
    sumP=sum(P,2);%每行元素之和
    for i=1:L
        P(i,:)=P(i,:)/sumP(i);
    end
    %M step
    MIUtip=MIU;
    SIGMAtip=SIGMA;
    WEIGHTtip=WEIGHT;
    for k=1:orders
        theta1=0;
        theta2=0;
        for i=1:L
            theta1=theta1+P(i,k);
            theta2=theta2+P(i,k)*MixedGaussSeq(i);
        end
        MIU(k)=theta2/theta1;
        WEIGHT(k)=theta1/L;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
        theta=0;
        for i=1:L
            theta=theta+P(i,k)*(MixedGaussSeq(i)-MIU(k))^2;                    
        end                                                                                                      
        SIGMA(k)=theta/theta1;
    end
    theta=[sum(abs(MIU-MIUtip)),sum(abs(SIGMA-SIGMAtip)),sum(abs(WEIGHT-WEIGHTtip))];%参数改变量小于某个阈值时不再迭代
    delta=max(theta);
end
   



