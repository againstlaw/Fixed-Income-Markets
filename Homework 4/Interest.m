clc
clear all


% load data
HW = readtable('Homework 4.xlsx');
Paths=HW{:,:};

N=length(Paths(:,1));   % number of samples
T=length(Paths(1,:));
t=0:(T-1);

% zero-coupon bond
D=zeros(T,1);
D(1)=1;

for i=2:T
    r=mean(Paths(:,1:i),2);
    D(i)=mean(exp(-r*t(i)));
end


figure(1)
plot(t,D,'-ko','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)
grid on
title('Discount Factor','Fontsize',16)
xlabel('Matuity (years)', 'Fontsize',14)
ylabel(('{\it D}(T)'), 'Fontsize',14)


% five-year interest rate cap
r5=mean(Paths,2);
K=0.045;
PayOff=max(0,Paths(:,end)-K);
Price_cap=mean(exp(-r5*t(end)).*PayOff)


% five-year interest rate floor
K2=0.067;
PayOff2=max(0,K2-Paths(:,end));
Price_floor=mean(exp(-r5*t(end)).*PayOff2)


% which one is more valuable: call(caplet) or put(floorlet)
PayOff3=max(0,Paths(:,end)-K2);
Price_cap3=mean(exp(-r5*t(end)).*PayOff3)


% which one is more valuable: 
% call on short interest rate or call on average short interest rate
K3=0.063;
PayOff4=max(0,Paths(:,end)-K3);
Price_cap4=mean(exp(-r5*t(end)).*PayOff4)

PayOff5=max(0,r5-K3);
Price_cap_avg=mean(exp(-r5*t(end)).*PayOff5)


% standard deviation
SD=std(Paths(:,end))
SD_avg=std(r5)


