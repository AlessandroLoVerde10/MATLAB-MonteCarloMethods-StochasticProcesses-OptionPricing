%  Control variates

 % calcolo del'europea asiatica
clc 
clear
close all

%Parametri del modello

S0 = 100
K = 100
vol = 0.2
r = 0.01

T = 1
m1  = 10^5

t = 360

dt = T/t

nsim = 10^2


tt = linspace(0, T, t);

% m1 osservazioni

S1 = zeros(m1,t);

S1(:,1) = S0 ;  %sol. esatta
C = zeros(1,nsim);
Q1 = zeros(m1,t);
Q1(:,1) = S0;
for j = 1: nsim
for i =1:t-1
   dW = sqrt(dt) * randn(1,m1);
    S1(:,i+ 1) = S1(:,i) .* exp( (r-0.5*vol^2)*(dt) + vol*dW');    
end

    Q_T1(:,j) = sqrt(mean(S1(:,1:end).^2,2));
    Q_T2(:,j) = geomean(S1(:,1:end)');
    Q_T3(:,j) = S1(:,end);
  %  payoffP1 = K-Q_T1;
    payoffC1 = Q_T1 - K;
     payoffC2 = Q_T2 - K;
      payoffC3 = Q_T3 - K;
  %  P1 = exp(-r*T)*mean(max(payoffP1,0));
    C = exp(-r*T)*(max(payoffC1,0));
    C1 = exp(-r*T)*mean(max(payoffC1,0));
     C2 = exp(-r*T)*mean(max(payoffC2,0));
     C3 = exp(-r*T)*mean(max(payoffC3,0));
end

%PMean1 = mean(P1)
Cmean1 = mean(C1)
Cmean2 = mean(C2)
Cmean3 = mean(C3)
%Pstdev1 = std(P1)
Cstdev1 = std(C1)
Cstdev = std(C(1,:))/10
%ErrP1 = P1 - PMean1
ErrC1 = C1 - Cmean1;

StartDates = '12-March-2014';
EndDates = '12-March-2020';
Rates = r;   
Compounding = -1;
Basis = 1;
Settle = '12-March-2014';
ExerciseDates = '12-March-2015';
StockSpec = stockspec(vol, S0);
RateSpec = intenvset('ValuationDate', StartDates, 'StartDates', StartDates, ...
   'EndDates', EndDates, 'Rates', Rates, 'Compounding', ...
    Compounding, 'Basis', Basis);
price_g = asiansensbykv(RateSpec, StockSpec, 'call', K, Settle, ExerciseDates)
[Call_E,Put] = blsprice(S0,K,r,T,vol)

% da aggiustare!!
MCov = cov(C1,C3)
MCov2 = cov(C1,C2)
beta = MCov(1,2)/MCov(2,2)
beta2 = MCov2(1,2)/MCov2(2,2)
%    cv_mean = (np.exp(-r*T[j])*SS[:,-1]).mean()
 %   C_MC[k,j] = np.mean(Payoff) - beta * (cv_mean - S0)
 
CALLCV = C1 + beta*(Call_E - C3) ;

CALLCVmean = mean (CALLCV)

CALLCVerror = CALLCV - CALLCVmean ;

CALLCVstdev = std (CALLCV)

CALLCV2 = C1 + beta*(price_g - C2) ;

CALLCVmean2 = mean (CALLCV2)

CALLCVerror2 = CALLCV2 - CALLCVmean2 ;

CALLCVstdev2 = std (CALLCV2)

figure
subplot (1,2,1)
plot(C1, 'linewidth', 2)
xlabel('Simlations')
ylabel('Price Call')
axis square

hold on
plot(CALLCV, 'linewidth', 2)
legend('MC', 'Control variates')


subplot (1,2,2)
plot(ErrC1, 'linewidth', 2)
xlabel('Simulations')
ylabel('error Call')
axis square

hold on
plot(CALLCVerror , 'linewidth', 2)
legend('MC', 'Control variates')

