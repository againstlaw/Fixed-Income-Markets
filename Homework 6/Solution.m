clc
clear all

% [num1, txt1] = xlsread('Homework 6 voldat.xlsx');
% [num2, txt2] = xlsread('Homework 6 pfilea.xlsx');

% vol=num1;
% D=num2;

% load data
filename1 = 'Homework 6 voldat.xlsx';
opts1 = detectImportOptions(filename1,'NumHeaderLines',0);
Volatility = readtable(filename1,opts1);

filename2 = 'Homework 6 pfilea.xlsx';
opts2 = detectImportOptions(filename2,'NumHeaderLines',0);
DT = readtable(filename2,opts2);

Volatility.Properties.VariableNames = {'Volatility'};
DT.Properties.VariableNames = {'DT'};

vol=Volatility{:,1};
D=DT{:,1};

N=length(D);

Matuity=0.5:0.5:N*0.5;


figure(1)
subplot(1,2,1)
plot(Matuity,D,'-ko','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)
grid on
title('Discount Factor','Fontsize',16)
xlabel('Matuity (years)', 'Fontsize',14)
ylabel('{\it D}(T)', 'Fontsize',14)

subplot(1,2,2)
plot(Matuity(2:N),vol,'-ko','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)
grid on
title('Volatility of Short Rate','Fontsize',16)
xlabel('Matuity (years)', 'Fontsize',14)
ylabel('\sigma', 'Fontsize',14)


% spot rate
SpotRate=zeros(N,1);

for i=1:N
    SpotRate(i)=2*(D(i)^(-1/2/Matuity(i))-1);
end


% annual rate
r=zeros(N,N);
r(1,1)=2*(1/D(1)-1);


for i=2:N
    dt=0.5;
    r_old=r(1:i-1,1:i-1);
    f = @(x) myFun(x,r_old,D(i),vol(i-1),i);
    x0=0.1;
    x = fsolve(f,x0);
    r(i,1)=x;   
    ratio=exp(-2*vol(i-1)*sqrt(dt));
    for j=2:i
        r(i,j)=r(i,j-1)*ratio;
    end
end

% expected value of r
prob=zeros(N,N);
for i=1:N
    ii=i-1;
    for j=1:i
        jj=j-1;
        prob(i,j)=factorial(ii)/factorial(jj)/factorial(ii-jj)*0.5^ii;
    end
end

ER=sum(r.*prob,2);

myD=zeros(N,1);
% for i=1:N
%     myD(i)=exp(sum(-ER(1:i)*0.5));
% end

myr=zeros(N,1);
for i=1:N
    myD(i)=1/prod(1+ER(1:i)*0.5);
    myr(i)=2*(myD(i)^(-1/2/Matuity(i))-1);
end

figure(2)
subplot(1,2,1)
plot(Matuity,D,'-ko','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)
hold on
plot(Matuity,myD,'b','linewidth',1)
grid on
title('Discount Factor','Fontsize',16)
xlabel('Matuity (years)', 'Fontsize',14)
ylabel('{\it D}(T)', 'Fontsize',14)
legend({'Data','BDT'}, 'FontSize', 12, 'Location','northeast')

subplot(1,2,2)
plot(Matuity,SpotRate,'-ko','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)
hold on
plot(Matuity,myr,'b','linewidth',1)
grid on
title('Spot Rate','Fontsize',16)
xlabel('Matuity (years)', 'Fontsize',14)
ylabel('{\it r}', 'Fontsize',14)
legend({'Data','BDT'}, 'FontSize', 12, 'Location','southeast')



function F = myFun(x,r,DT,vol,N)

% semi-annual
dt=0.5;
ratio=exp(-2*vol*sqrt(dt));

% short rate at step N
myR=zeros(N,N);
myR(1:N-1,1:N-1)=r;   % spot rate at previous steps
myR(N,1)=x;
for i=2:N
    myR(N,i)=myR(N,i-1)*ratio;
end

P=zeros(N+1,N+1);
P(N+1,:)=1;
for n=N:-1:1
    for m=1:n
        P(n,m)=(0.5*P(n+1,m)+0.5*P(n+1,m+1))/(1+myR(n,m)/2);
    end
end

F = DT-P(1,1);

end


