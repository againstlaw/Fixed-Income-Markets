clc
clear all


% load data
Data = readtable('Homework_5.xlsx','Sheet','Sheet1');

DateMatrix=Data{:,1:3};
myDate=arrayfun(@datestr,datenum(DateMatrix),'UniformOutput',false);
Date=datetime(myDate);

CMTList={'CMT(0.25)','CMT(2)','CMT(3)','CMT(5)','CMT(7)','CMT(10)'};


x0=[0.005 0.5 0.01 0.1 0.02];

InputData=Data{:,5:10}/100;

f = @(x) myFun(x,InputData);
% [x,fval] = fminunc(f,x0)

options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton','MaxFunctionEvaluations',2000);
[x,fval,exitflag,output] = fminunc(f,x0,options)


[X,Y]=FunXY(x,InputData);

x(1)/x(2)
EstimationMean=mean(X)

x(3)^2/2/x(2)
EstimationSD=std(X)


r=X+Y;


figure(1)
plot(Date,Data{:,5:10},'-','linewidth',1)
datetick('x')
grid on
title('CMT','Fontsize',14)
xlabel('Date', 'Fontsize',14)
ylabel('YTM (%)', 'Fontsize',14)
legend(CMTList, 'FontSize', 12, 'Location','north','NumColumns',2)


figure(2)
plot(Date,r,'-r','linewidth',1)
datetick('x')
grid on
xlabel('Date', 'Fontsize',14)
title('r=X+Y','Fontsize',14)




YTM_computed = FunYTM(x,InputData);

error=InputData(:,2:5)-YTM_computed;


figure(3)
subplot(2,1,1)
plot(Date,error(:,1),'-r','linewidth',1)
datetick('x')
grid on
xlabel('Date', 'Fontsize',14)
title('CMT(2)','Fontsize',14)

subplot(2,1,2)
plot(Date,error(:,2),'-b','linewidth',1)
datetick('x')
grid on
xlabel('Date', 'Fontsize',14)
title('CMT(3)','Fontsize',14)



figure(4)
subplot(2,1,1)
plot(Date,error(:,3),'-r','linewidth',1)
datetick('x')
grid on
xlabel('Date', 'Fontsize',14)
title('CMT(5)','Fontsize',14)

subplot(2,1,2)
plot(Date,error(:,4),'-b','linewidth',1)
datetick('x')
grid on
xlabel('Date', 'Fontsize',14)
title('CMT(7)','Fontsize',14)



% Time Series Properties
% AR(3)
mySample=error(:,3);
Model = arima(3,0,0);
EstMdl = estimate(Model,mySample);


figure(5)
subplot(1,2,1);
autocorr(mySample)
subplot(1,2,2)
parcorr(mySample)

% Plot the residuals
res = infer(EstMdl,mySample);
stdRes = res/sqrt(EstMdl.Variance);    % Standardized residuals


figure(6)
plot(stdRes);
title 'Inferred Residuals';
xlabel('numPeriods', 'Fontsize',10)


figure(7)
autocorr(stdRes)

[h,pValue] = lbqtest(stdRes,'lags',[1,2])






function RMSE = myFun(x,data)

alphaX=x(1);
betaX=x(2);
sigmaX=x(3);
alphaY=0;
betaY=x(4);
sigmaY=x(5);

N=length(data(:,1));
T1=0.25;
T2=10;
AX1=A(alphaX,betaX,sigmaX,T1);
BX1=B(betaX,T1);
AY1=A(alphaY,betaY,sigmaY,T1);
BY1=B(betaY,T1);
AX2=A(alphaX,betaX,sigmaX,T2);
BX2=B(betaX,T2);
AY2=A(alphaY,betaY,sigmaY,T2);
BY2=B(betaY,T2);

a11=BX1/T1;
a12=BY1/T1;
a21=BX2/T2;
a22=BY2/T2;
AA=[a11 a12; a21 a22];

for i=1:N
    b1=data(i,1)+log(AX1)/T1-log(AY1)/T1;
    b2=data(i,6)+log(AX2)/T2-log(AY2)/T2;
    b=[b1; b2];
    result=AA\b;
    X(i,1)=result(1);
    Y(i,1)=result(2);
end


Maturity=[2 3 5 7];
M=length(Maturity);

% compute Discount rate
myT=0.5:0.5:10;
DT_computed=DT(alphaX,betaX,sigmaX,alphaY,betaY,sigmaY,X,Y,myT);
YTM_computed=zeros(N,M);

% compute par rate
for j=1:M
    T=Maturity(j);
    num=2*T;
    ParRate=2*(1-1*DT_computed(:,num))./sum(DT_computed(:,1:num),2);
    YTM_computed(:,j)=ParRate;
end

YTM=data(:,2:5);

% compute the root-mean-squared-error
err=abs(YTM_computed-YTM);
RMSE=sqrt(sum(sum(err.^2)))/4;

end




function [X,Y] = FunXY(x,data)

alphaX=x(1);
betaX=x(2);
sigmaX=x(3);
alphaY=0;
betaY=x(4);
sigmaY=x(5);

N=length(data(:,1));
T1=0.25;
T2=10;
AX1=A(alphaX,betaX,sigmaX,T1);
BX1=B(betaX,T1);
AY1=A(alphaY,betaY,sigmaY,T1);
BY1=B(betaY,T1);
AX2=A(alphaX,betaX,sigmaX,T2);
BX2=B(betaX,T2);
AY2=A(alphaY,betaY,sigmaY,T2);
BY2=B(betaY,T2);

a11=BX1/T1;
a12=BY1/T1;
a21=BX2/T2;
a22=BY2/T2;
AA=[a11 a12; a21 a22];

for i=1:N
    b1=data(i,1)+log(AX1)/T1-log(AY1)/T1;
    b2=data(i,6)+log(AX2)/T2-log(AY2)/T2;
    b=[b1; b2];
    result=AA\b;
    X(i,1)=result(1);
    Y(i,1)=result(2);
end


end



function myDT = DT(alphaX,betaX,sigmaX,alphaY,betaY,sigmaY,X,Y,T)

AX=A(alphaX,betaX,sigmaX,T);
BX=B(betaX,T);
AY=A(alphaY,betaY,sigmaY,T);
BY=B(betaY,T);

myDT = AX.*AY.*exp(-BX.*X-BY.*Y);

end



function YTM_computed = FunYTM(x,data)

alphaX=x(1);
betaX=x(2);
sigmaX=x(3);
alphaY=0;
betaY=x(4);
sigmaY=x(5);

N=length(data(:,1));
T1=0.25;
T2=10;
AX1=A(alphaX,betaX,sigmaX,T1);
BX1=B(betaX,T1);
AY1=A(alphaY,betaY,sigmaY,T1);
BY1=B(betaY,T1);
AX2=A(alphaX,betaX,sigmaX,T2);
BX2=B(betaX,T2);
AY2=A(alphaY,betaY,sigmaY,T2);
BY2=B(betaY,T2);

a11=BX1/T1;
a12=BY1/T1;
a21=BX2/T2;
a22=BY2/T2;
AA=[a11 a12; a21 a22];

for i=1:N
    b1=data(i,1)+log(AX1)/T1-log(AY1)/T1;
    b2=data(i,6)+log(AX2)/T2-log(AY2)/T2;
    b=[b1; b2];
    result=AA\b;
    X(i,1)=result(1);
    Y(i,1)=result(2);
end


Maturity=[2 3 5 7];
M=length(Maturity);

% compute Discount rate
myT=0.5:0.5:10;
DT_computed=DT(alphaX,betaX,sigmaX,alphaY,betaY,sigmaY,X,Y,myT);
YTM_computed=zeros(N,M);

% compute par rate
for j=1:M
    T=Maturity(j);
    num=2*T;
    ParRate=2*(1-1*DT_computed(:,num))./sum(DT_computed(:,1:num),2);
    YTM_computed(:,j)=ParRate;
end


end



function Yield = CMT(alphaX,betaX,sigmaX,alphaY,betaY,sigmaY,X,Y,T)

AX=A(alphaX,betaX,sigmaX,T);
BX=B(betaX,T);
AY=A(alphaY,betaY,sigmaY,T);
BY=B(betaY,T);

Yield = -log(AX)/T+BX/T*X+log(AY)/T+BY/T*Y;

end



function myA = A(alpha,beta,sigma,T)

myA=(sigma^2/2/beta^2-alpha/beta)*T+(alpha/beta^2-sigma^2/beta^3)*(1-exp(-beta*T))...
    +sigma^2/4/beta^3*(1-exp(-2*beta*T));
myA=exp(myA);

end



function myB = B(beta,T)

myB=1/beta*(1-exp(-beta*T));

end


