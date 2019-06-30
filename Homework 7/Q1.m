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

Matuity = 0.5:0.5:N*0.5;


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


L = chol(rho,'lower');


% Forward Par rates for semiannual coupon bonds
FT=[1 2 3 4 5];
M=length(FT);
NT=5*2;
ParRate=zeros(M,1);

for i=1:M
    pos=FT(i)*2;
    ParRate(i)=2*(100*D(NT)-100*D(NT+pos))/sum(D(NT+1:NT+pos));
end


figure(2)
plot(FT,ParRate,'-ko','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)
grid on
title('Forward Par Curve','Fontsize',16)
xlabel('Matuity (years)', 'Fontsize',14)
ylabel(('{\it c} (%)'), 'Fontsize',14)


% T = 20 periods
T = length(L(:,1));
Nsim = 100000;

vol = [0; vol];   % sigma(0.5)=0

myD = zeros(T,T,Nsim);
dt = 0.5;

for k=1:Nsim
    RandomNum = L*randn(T,T);
    sample = zeros(T,T);
    old_D=D;
    for j=1:T
        r = 2*(1/old_D(j)-1);
        new_D = zeros(T,1);
        for i=j:T
            new_D(i) = old_D(i)+r*old_D(i)*dt + vol(i-j+1)*RandomNum(i,j)*sqrt(dt);
        end
        sample(:,j)=new_D;
        old_D=new_D;
    end
    myD(:,:,k)=sample;
    
end

myT=0.5:0.5:10;

% sample spot rate
SpotRate=zeros(T,T);
for j=1:T
    for i=j:T
        SpotRate(i,j)=2*(sample(i,j)^(-1/2/myT(i))-1);
    end
end


DD=mean(myD,3);

% t = 5 years
DDDDD=DD(NT+1:2*NT,NT);
P=zeros(M,1);

for i=1:M
    pos=FT(i)*2;
    P(i)=ParRate(i)/2*sum(DDDDD(1:pos))+100*DDDDD(pos);
end
    

FiveYear=zeros(NT,NT);
for i=1:NT
    FiveYear(:,i)=DD(i+1:i+NT,i);
end


% 3-D Plot
newT=0.5:0.5:5;
figure(3)
[X,Y]= meshgrid(newT,newT);
h=surf(X,Y,FiveYear);
% set(h,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
% surf(Date,T,myData','EdgeColor','none','LineStyle','none','FaceLighting','phong')
colorbar
xlabel('Date (Year)', 'Fontsize',14)
ylabel('Maturity (Year)', 'Fontsize',14)
zlabel('{\it D}(T)', 'Fontsize',14)
title('Sample Average','Fontsize',16)

