clc
clear all


% load data
HW = readtable('Homework 2 Data.xlsx','Sheet','T-Note');
Price=HW{:,1};
Matuity=HW{:,2};
Coupon=HW{:,3};
Yield=HW{:,4};


Y=Yield;
T=Matuity;
T2=T.*T;
T3=T.*T.*T;
T4=T.*T.*T.*T;
T5=T.*T.*T.*T.*T;

X=[ones(size(T)) T T2 T3 T4 T5];

b=regress(Y,X);

myT=0.5:0.5:25;
myY=b(1)+b(2)*myT+b(3)*myT.^2+b(4)*myT.^3+b(5)*myT.^4+b(6)*myT.^5;


figure(1)
plot(Matuity,Y,'o','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)
hold on
plot(myT,myY,'r','linewidth',1)
grid on
title('Smoothing Yield Curve','Fontsize',16)
xlabel('Matuity (years)', 'Fontsize',14)
ylabel('{\it Y}(T) (%)', 'Fontsize',14)
legend({'T-Note data','Regression'}, 'FontSize', 12, 'Location','southeast')


% Bootstraping
% Semiannual-pay bond
N=length(myT);
FaceValue=100;
D=zeros(N,1);
D(1)=FaceValue/(0.5*myY(1)+FaceValue);   % discount factor

for i=2:N
    coupon=0.5*myY(i)/100*FaceValue;
    D(i)=(FaceValue-coupon*sum(D(1:i-1)))/(coupon+FaceValue);
end

r=zeros(N,1);   % spot rate
for i=1:N
    r(i)=2*(D(i)^(-1/2/myT(i))-1);
end


% 6-month Forward Rates
m=6/12;
ForwardRate=zeros(N-1,1);

for i=1:N-1
    ForwardRate(i)=1/m*(D(i)/D(i+1)-1);
end

figure(2)
plot(myT,r,'r','linewidth',1)
hold on
plot(myT,myY/100,'b','linewidth',1)
hold on
plot(myT(1:N-1),ForwardRate,'k','linewidth',1)
grid on
title('Bootstrapping Curves','Fontsize',16)
xlabel('Matuity (years)', 'Fontsize',14)
legend({'Spot Curve','Par Curve', '6-month Forward'}, 'FontSize', 12, 'Location','southeast')    




