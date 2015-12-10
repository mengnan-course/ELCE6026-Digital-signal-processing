clear all;
close all;
addpath('ARMAX_GARCH_K_Toolbox');
% Reset random stream for reproducibility.
rng(0,'twister');
JAssets = 10; % desired number of assets
% Generate means of returns between -0.02 and 0.05.
a = -0.01; b = 0.05;
mean_return = a + (b-a).*rand(JAssets,1);
% Generate standard deviations of returns between 0.008 and 0.006.
a = 0.08; b = 0.06;
stdDev_return = a + (b-a).*rand(JAssets,1);
%%
Ntime=200;
%% X: Each row is a time-instant. Each Column is an asset.
X=zeros(Ntime,JAssets);
for j=1:JAssets
    X(:,j)=mean_return(j)+stdDev_return(j)*randn(Ntime,1);
end
%% 
%% MF: mean forecast
%% VF: forecast of variance
MFpred=zeros(JAssets,1);
VFpred=zeros(JAssets);
for j=1:JAssets
    %% For each of the variables, fit the ARMA(1,1)-GARCH(1,1) model
    [parameters, stderrors, LLF, ht, resids, summary] = garch(X(:,j),'GARCH', 'GAUSSIAN',1,1,0,1,1,0,[]);
    %% 1-step ahead Prediction of the mean and covariance of return 
    [MFpred(j), VFpred(j,j), ~, ~] = garchfor2(X(:,j), resids, ht, parameters, 'GARCH', 'GAUSSIAN',1,1,1,1,1);
end


a=ones(JAssets,1)'*(VFpred\ones(JAssets,1));
b=ones(JAssets,1)'*(VFpred\MFpred);
c=MFpred'*(VFpred\MFpred);

target=[-.1:5e-3:.1]';
risk=zeros(length(target),1);
for j=1:length(target)
   delta=a*c-b^2;
    lambda1=(c-b*target(j))/delta;
    lambda2=(a*target(j)-b)/delta;
    w=VFpred\(lambda1*ones(JAssets,1)+lambda2*MFpred);
    risk(j)=w'*VFpred*w;
end
plot(risk,target);
ylabel('Expected return')
xlabel('Variance')

%% Write your own code;
Delta = 4*10^-5;
delta1 = 10^-3;
index1 = find(abs(risk - delta1) < Delta);
mu_1 = max(target(index1));

delta2 = 1.6 * 10^-3;
index2 = find(abs(risk - delta2) < Delta);
mu_2 = max(target(index2));

hold on
plot([delta1,delta2],[mu_1,mu_2],'r*')
text(delta1,mu_1,['\delta_1 = ' num2str(delta1) ' \mu_1 = ' num2str(mu_1)])
text(delta2,mu_2,['\delta_2 = ' num2str(delta2) ' \mu_2 = ' num2str(mu_2)])


%% find the better one 
index_1 = find(target == mu_1)
index_2 = find(target == mu_2)

Y1 = mu_1;
X1 = risk(index_1);
Y2 = mu_2;
X2 = risk(index_2);

Y1_1 = target(index_1 + 1);
X1_1 = risk(index_1 + 1);
Y2_1 = target(index_2 + 1);
X2_1 = risk(index_2 + 1);


slope_1 = ((Y1_1+Y1)/2 - Y1) / ((X1_1+X1)/2 - X1);
slope_2 = ((Y2_1+Y2)/2 - Y2) / ((X2_1+X2)/2 - X2);

RangeX = max(risk) - min(risk);
X1_start = min(risk);
Y1_start = Y1 + slope_1*(X1_start - X1);
X1_end = max(risk);
Y1_end = Y1 + slope_1*(X1_end - X1);

X2_start = min(risk);
Y2_start = Y2 + slope_2*(X2_start - X2);
X2_end = max(risk);
Y2_end = Y2 + slope_2*(X2_end - X2);

a = plot([X1_start,X1_end],[Y1_start,Y1_end],'g-','LineWidth',1);
b = plot([X2_start,X2_end],[Y2_start,Y2_end],'r-','LineWidth',1);
axis([0.4*10^-3, 2.2*10^-3, -0.1, 0.1])

L1 = ['slope1 = ' num2str(slope_1)];
L2 = ['slope2 = ' num2str(slope_2)];
legend([a,b],L1,L2,'Location','southeast');






