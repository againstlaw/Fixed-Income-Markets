clc
clear all


% load data
HW = readtable('Homework 2 Data.xlsx');
Matuity=HW{:,1};
Price=HW{:,2};

FaceValue=100;
D=Price/FaceValue;   % discount factor

N=length(Matuity);

% spot rate
r=zeros(N,1);

for i=1:N
    r(i)=2*(D(i)^(-1/2/Matuity(i))-1);
end

% 3-month forward rate
m=3/12;
f=zeros(N-1,1);

for i=1:N-1
    f(i)=1/m*(D(i)/D(i+1)-1);
end


figure(1)
subplot(1,2,1)
plot(Matuity,r,'linewidth',1)
grid on
title('Spot Curve','Fontsize',16)
xlabel('Matuity (years)', 'Fontsize',14)

subplot(1,2,2)
plot(Matuity(1:N-1),f,'linewidth',1)
grid on
title('3-month Forward Curve','Fontsize',16)
xlabel('Matuity (years)', 'Fontsize',14)


T=Matuity;
T2=T.*T;
T3=T.*T.*T;
T4=T.*T.*T.*T;
T5=T.*T.*T.*T.*T;

X=[T T2 T3 T4 T5];
Y=log(D);
b=regress(Y,X);

DD=exp(b(1)*T+b(2)*T2+b(3)*T3+b(4)*T4+b(5)*T5);

figure(2)
plot(Matuity,D,'-ko','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)
hold on
plot(Matuity,DD,'r','linewidth',1)
grid on
title('Spot Curve','Fontsize',16)
xlabel('Matuity (years)', 'Fontsize',14)
ylabel('{\it D} (T)', 'Fontsize',14)
legend({'STRIPS data','Regression'}, 'FontSize', 12, 'Location','southeast')


% Problem 3
myT=0.5:0.5:25;
myD=exp(b(1)*myT+b(2)*myT.^2+b(3)*myT.^3+b(4)*myT.^4+b(5)*myT.^5);
M=length(myT);

% spot rate
myr=zeros(M,1);

for i=1:M
    myr(i)=2*(myD(i)^(-1/2/myT(i))-1);
end

figure(3)
plot(myT,myr,'-ko','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)
grid on
title('Spot Curve','Fontsize',16)
xlabel('Matuity (years)', 'Fontsize',14)



% Problem 4
% Par rate
ParRate=zeros(M,1);

for i=1:M
    ParRate(i)=2*(100-100*myD(i))/sum(myD(1:i));
end

figure(4)
plot(myT,ParRate,'-ko','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)
grid on
title('Par Curve','Fontsize',16)
xlabel('Matuity (years)', 'Fontsize',14)
ylabel(('{\it c} (%)'), 'Fontsize',14)


% Problem 5
% 6-month Forward Rates
m=6/12;
newT=0.5:0.5:25.5;
newD=exp(b(1)*newT+b(2)*newT.^2+b(3)*newT.^3+b(4)*newT.^4+b(5)*newT.^5);

ForwardRate=zeros(M,1);

for i=1:M
    ForwardRate(i)=1/m*(newD(i)/newD(i+1)-1);
end


figure(5)
plot(myT,ForwardRate,'-ko','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)
grid on
title('6-month Forward Curve','Fontsize',16)
xlabel('Matuity (years)', 'Fontsize',14)

