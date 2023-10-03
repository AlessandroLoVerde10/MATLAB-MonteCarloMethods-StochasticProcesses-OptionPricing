% SABR simulation
 % calcolo del'europea asiatica
%% Real call prices
clear all
clc
close all
[num] = xlsread('Implied volatility call FTSEMIB 19MAY21')
r = -0.05
S0 = 24486.45;
MarketStrikes = num(2:end,1)'

MarketVolatilit = num(2:end, 2:end)/100
MarketVolatilit(6,7) = 0.178
MarketVolatilit(7,6) = 0.175
MarketVolatilities = MarketVolatilit
%CurrentForwardValues = [24377.57; 24221.71; 23929.08; 23807.87; 23344.49; 22971.68; 22417.75; 22079.33]
CurrentForwardValues = [S0; S0; S0; S0; S0; S0; S0; S0]
NumMaturities = 8
Beta = 1

Settle = datenum('19-May-2021');

%ExerciseDate = [Settle + (364/12) ; Settle + (364/12)+ (364/4); Settle + (364/12)+ (364/2); Settle + (364/12)+(364/4)*3; Settle + (364/12)+ 364; Settle + (364/12)+ (364/4)*6; Settle + (364/12)+ (364/4)*8; Settle + (364/12)+ (364/4)*12]
ExerciseDate = [datenum('18-Jun-2021') ; datenum('17-Sep-2021'); datenum('17-Dec-2021'); datenum('18-Mar-2022'); datenum('17-Jun-2022'); datenum('16-Dec-2022'); datenum('16-Jun-2023'); datenum('15-Dec-2023')]
YearsToExercise = yearfrac(Settle, ExerciseDate, 1);
for i = 1:8
CurrentForwardValues1 (i) = exp(r*YearsToExercise(i,1))*S0
end
Call_E = zeros(8,8)
for j = 1:8
for i = 1:8
[Call_E(i, j),PUT] = blsprice(S0, MarketStrikes(i) ,r,YearsToExercise (j),MarketVolatilities(i,j))
end
end
MarketStrike = repmat(MarketStrikes, 8, 1)'
Betas = repmat(Beta, NumMaturities, 1)
Alphas = zeros(NumMaturities, 1);
Rhos = zeros(NumMaturities, 1);
Nus = zeros(NumMaturities, 1);

%options = optimoptions('lsqnonlin','Display','none');

for k = 1:NumMaturities
    % Fit Rho and Nu (while converting at-the-money volatility into Alpha)
    objFun = @(X) MarketVolatilities(:,k) -  blackvolbysabr(X(1), Betas(k),  X(2),  X(3), Settle, ExerciseDate(k), CurrentForwardValues(k), MarketStrike(:,k));
    
    X = lsqnonlin(objFun, [0.5 0 0.5], [0 -1 0], [Inf 1 Inf]);
    
    Alphas(k) = X(1);
    Rhos(k) = X(2);
    Nus(k) = X(3);
end

CalibratedPrameters = [Alphas Betas Rhos Nus]
%%
%Parametri del modello
%0.2100    1.0000   -0.7279    1.3779
a=5
r = -0.05;
S0 = 24486.45;
%    1.0000      
alpha0 =  CalibratedPrameters(a,1)  ;
beta = CalibratedPrameters(a,2);
rho = CalibratedPrameters(a,3) ;
nu =  CalibratedPrameters(a,4);

T = YearsToExercise(a)
m2  = 10^3
barriera = 21000;
t = 400
dt = T/t;
nsim = 100;
K = [21000, 22000, 23000, 24000, 25000, 26000, 27000, 28000];
tt = linspace(0, T, t);
% m1 osservazioni
SS2= zeros(m2,t);
SS= zeros(m2,t);
alphas = zeros(m2,t);
SS2(:,1) = S0 ;  %sol. esatta
alphas(:,1) = alpha0;
SS(:,1) = S0 ; 

for j = 1: nsim
 for i = 1:t-1
   dW1 = sqrt(dt) * randn(1,m2);
   dW2 = sqrt(dt) * randn(1,m2);
   alphas(:,i+1) = alphas(:,i) + nu*alphas(:,i).*dW2';
   SS2(:,i+1) = SS2(:,i) + r*SS2(:,i)*dt + (alphas(:,i).*(SS2(:,i).^beta).*(sqrt(1-rho^2)*dW1+rho*dW2)');
   for p = 1:m2-1
   if SS(p,i) <= barriera
         SS(p, i+1) = 0;
   else
     SS(p,i+1) = SS(p,i) + r*SS(p,i)*dt + (alphas(p,i).*(SS(p,i).^beta).*(sqrt(1-rho^2)*dW1(p)+rho*dW2(p))');
   end
   end 

 end
 
   for k = 1:numel(K)

  Q_T2(:,j) = SS2(:,end);
  Q_T(:,j) = SS(:,end);
 payoffP2 = K(k)-Q_T2;
 payoffC2 = Q_T2 - K(k);
  payoffP = K(k)-Q_T;
 payoffC = Q_T - K(k);
 Cx = exp(-r*T)*(max(payoffC2(:,1),0));
 CS = std(Cx)/sqrt(m2);
 P2 = exp(-r*T)*mean(max(payoffP2,0));
 C2 = exp(-r*T)*mean(max(payoffC2,0));
 P = exp(-r*T)*mean(max(payoffP,0));
 C = exp(-r*T)*mean(max(payoffC,0));

Pmean2(k) = mean(P2);
Cmean2(k) = mean(C2);

Pstdev2 (k)= std(P2);
Cstdev2 (k) = std(C2);
%ErrP2 (:,k) = P2 - Pmean2(k);
%ErrC2(:,k) = C2 - Cmean2(k);

Pmean (k)= mean(P);
Cmean(k) = mean(C);
 end
end
%figure

%plot(SS', 'linewidth', 1)
%xlabel('Time')
%ylabel('Price S')
%axis square

%hold on

%yline(barriera,'-',{'Barrier'}, 'linewidth', 2)



%subplot (2,2,2)
%plot(P1, 'linewidth', 2)
%xlabel('Time')
%ylabel('Price Put')
%axis square

%subplot (2,2,3)
%plot(ErrC1, 'linewidth', 2)
%xlabel('Time')
%ylabel('Price Call')
%axis square

%subplot (2,2,4)
%plot(ErrP1, 'linewidth', 2)
%xlabel('Time')
%ylabel('Price Put')
%axis square

figure 
d = plot (K, Cmean, 'b', linewidth = 2)
d(1).Marker = '*';

hold on 
f = plot (K, Cmean2, 'r', linewidth = 2)
f(1).Marker = '*';
hold on
p = plot (K, Call_E(1:8,a), 'g', linewidth = 2)
p(1).Marker = '*';


xlabel('Strikes', 'Fontweight', 'bold');
ylabel('Call Prices', 'Fontweight', 'bold');
legend('Barrier Call prices simulated with SABR', 'EU Call prices simulated with SABR','Observed Call prices')


