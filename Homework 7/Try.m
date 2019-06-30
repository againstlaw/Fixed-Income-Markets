clc
clear all


% load data
filename1 = 'Homework 7 sigma.xlsx';
opts1 = detectImportOptions(filename1,'NumHeaderLines',0);
Volatility = readtable(filename1,opts1);

filename2 = 'Homework 7 pfilea.xlsx';
opts2 = detectImportOptions(filename2,'NumHeaderLines',0);
DT = readtable(filename2,opts2);

filename3 = 'Homework 7 corrin.csv';
opts3 = detectImportOptions(filename3,'NumHeaderLines',0);
Correlation = readtable(filename3,opts3);

Volatility.Properties.VariableNames = {'Volatility'};
DT.Properties.VariableNames = {'DT'};

vol = Volatility{:,1};
D = DT{:,1};
rho = Correlation{:,:};

N = length(D);

Maturity = 0.5:0.5:N*0.5;


figure(1)
subplot(1,2,1)
plot(Maturity,D,'-ko','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)
grid on
title('Discount Factor','Fontsize',16)
xlabel('Maturity (years)', 'Fontsize',14)
ylabel('{\it D}(T)', 'Fontsize',14)

subplot(1,2,2)
plot(Maturity(2:N),vol,'-ko','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)
grid on
title('Volatility of Short Rate','Fontsize',16)
xlabel('Maturity (years)', 'Fontsize',14)
ylabel('\sigma', 'Fontsize',14)


L = chol(rho,'lower');



% T = 20 periods
T = length(L(:,1));
Nsim = 100000;

vol = [0; vol];   % sigma(0.5)=0

myD = zeros(T,T,Nsim);
dt = 0.5;
myr = zeros(T,Nsim);
r_sample=zeros(T,1);

for k=1:Nsim
    sample = zeros(T,T);
    old_D=D(1:T);
    for j=1:T
        r = 2*(1/old_D(j)-1);
        new_D = zeros(T,1);
        %RandomNum = L(j:T,j).*randn(T-j+1,1);
        RandomNum = L(1:T-j+1,j).*randn(T-j+1,1);
        for i=j:T
            new_D(i) = old_D(i) + r*old_D(i)*dt + vol(i-j+1)*RandomNum(i-j+1)*sqrt(dt);
        end
        sample(:,j)=new_D;
        old_D=new_D;
        r_sample(j)=r;
    end
    myD(:,:,k)=sample;
    myr(:,k)=r_sample;
end

DD=mean(myD,3);
rr=mean(myr,2);


figure(2)
plot(Maturity(1:20),rr,'-ko','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)
grid on
title('Average Interest Rate','Fontsize',16)
xlabel('Maturity (years)', 'Fontsize',14)
ylabel(('{\it r}'), 'Fontsize',14)

