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



% Problem 1
% Par rate for semiannual-coupon bonds
T=1:25;
M=length(T);
ParRate=zeros(M,1);

for i=1:M
    ParRate(i)=2*(100-100*D(2*i))/sum(D(1:2*i));
end


figure(2)
plot(T,ParRate,'-ko','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)
grid on
title('Par Curve','Fontsize',16)
xlabel('Matuity (years)', 'Fontsize',14)
ylabel(('{\it c} (%)'), 'Fontsize',14)



% Problem 2
% DV01
c=ParRate;
newYield=ParRate+0.01;
newYield2=ParRate-0.01;
P=zeros(M,1);
PP=zeros(M,1);

for i=1:M
    YTM=newYield(i)/100/2;
    P(i)=c(i)/2/YTM*(1-1/(1+YTM)^(2*i))+100/(1+YTM)^(2*i);
    YTM2=newYield2(i)/100/2;
    PP(i)=c(i)/2/YTM2*(1-1/(1+YTM2)^(2*i))+100/(1+YTM2)^(2*i);
end

DV01_up=100-P;
DV01_down=PP-100;


figure(3)
plot(T,DV01_up,'-ko','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)
grid on
title('DV01','Fontsize',16)
xlabel('Matuity (years)', 'Fontsize',14)



% Problem 3
% Macauley and modified durations
TT=1:5;
Macauley=zeros(5,1);
Modified=zeros(5,1);

for i=1:5
    CF=[c(i)/2*ones(1,2*i-1), c(i)+100];
    Years=cumsum(0.5*ones(1,2*i));
    Macauley(i)=Years.*CF*D(1:2*i)/100;
    halfyield=c(i)/2/100;
    Modified(i)=Macauley(i)/(1+halfyield);
end


figure(4)
plot(TT,Macauley,'-ko','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)
hold on
plot(TT,Modified,'-bs','linewidth',1,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',3)
grid on
title('Macauley Durations','Fontsize',16)
xlabel('Matuity (years)', 'Fontsize',14)
legend({'Macauley Durations','Modified Durations'}, 'FontSize', 12, 'Location','southeast')



% Problem 5
% convexity
Convexity=zeros(5,1);

for i=1:5
    CF=[c(i)/2*ones(1,2*i-1), c(i)+100];
    Periods=2*i;
    sum=0;
    for j=1:Periods
        sum=sum+j*(j+1)*CF(j)*D(j);
    end
    halfyield=c(i)/2/100;
    Convexity(i)=sum/(1+halfyield)^2/2^2/100;
end


figure(5)
plot(TT,Convexity,'-ko','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)
grid on
title('Convexity','Fontsize',16)
xlabel('Matuity (years)', 'Fontsize',14)


% Problem 6
% Price change in percentage
dy=100/100/100;   % 100 basis point

PriceChange=-Modified*dy*100+0.5*Convexity*dy^2*100;

% shited spot curve
Yield_shiftedup=ParRate+dy*100;
Yield_shifteddown=ParRate-dy*100;

Price_shiftedup=zeros(5,1);
Price_shifteddown=zeros(5,1);

for i=1:5
    YTM=Yield_shiftedup(i)/100/2;
    Price_shiftedup(i)=c(i)/2/YTM*(1-1/(1+YTM)^(2*i))+100/(1+YTM)^(2*i);
    YTM2=Yield_shifteddown(i)/100/2;
    Price_shifteddown(i)=c(i)/2/YTM2*(1-1/(1+YTM2)^(2*i))+100/(1+YTM2)^(2*i);
end

% price change in percentage
PriceChange_up=Price_shiftedup-100;
PriceChange_down=Price_shifteddown-100;


compare=[PriceChange,PriceChange_up];

figure(6)
bar(TT,compare)
xlabel('Matuity (years)', 'Fontsize',14)
ylabel('Price change (%)', 'Fontsize',14)
title('Upward Parallel Shift','Fontsize',16)
legend({'Approx','Actual'}, 'FontSize', 12, 'Location','southwest')
