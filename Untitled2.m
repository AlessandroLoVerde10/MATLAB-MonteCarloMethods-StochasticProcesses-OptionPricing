 % calcolo del'europea asiatica
clc 
clear
close all
%Parametri del modello
vol = 0.2
r = 0.01
S0 = 100
T = 1

m3 = 2*10^6

t = 5

dt = 1/t

nsim = 10^2
K = 100

tt = linspace(0, T, t);

% 10^ osservazioni

S3 = zeros(m3,t);

S3(:,1) = S0 ;  %sol. esatta

Q3 = zeros(m3,nsim);
Q3 (:,1) = S0;
for j = 1: nsim
for i =1:t-1
   dW = sqrt(dt) * randn(1,m3);
    S3(:,i+ 1) = S3(:,i) .* exp( (r-0.5*vol^2)*(dt) + vol*dW');
end
    Q_T3(:,j) = sqrt(mean(S3(:,1:end).^2,2));
  %  payoffP3 = K-Q_T3;
    payoffC3 = Q_T3 - K;
   % P3 = exp(-r*T)*mean(max(payoffP3,0));
    C3 = exp(-r*T)*mean(max(payoffC3,0));
end

%Pmean3 = mean(P3)
Cmean3 = mean(C3)
%Pstdev3 = std(P3)
Cstdev3 = std(C3)

%ErrP3 = P3 - Pmean3;
ErrC3 = C3 - Cmean3;


%% Antitetiche


m1  = 10^6;

S1 = zeros(m1,t);

S1(:,1) = S0 ;  %sol. esatta

S3 = zeros(m1,t);

S3(:,1) = S0 ; 

Q1 = zeros(m1,nsim);
Q1(:,1) = S0;
for j = 1: nsim
for i =1:t-1
   dW = sqrt(dt) * randn(1,m1);
    S1(:,i+ 1) = S1(:,i) .* exp( (r-0.5*vol^2)*(dt) + vol*dW');
      S3(:,i+ 1) = S3(:,i) .* exp( (r-0.5*vol^2)*(dt) - vol*dW');
    S = [S3;S1];
end
    Q_T1(:,j) = sqrt(mean(S(:,1:end).^2,2));
%    payoffP1 = K-Q_T1;
    payoffC1 = Q_T1 - K;
  %  P1 = exp(-r*T)*mean(max(payoffP1,0));
    C1 = exp(-r*T)*mean(max(payoffC1,0));
end

%PMean1 = mean(P1)
Cmean1 = mean(C1)
%Pstdev1 = std(P1)
Cstdev1 = std(C1)
%ErrP1 = P1 - PMean1
ErrC1 = C1 - Cmean1

figure
subplot (1,2,1)
plot(C1, 'linewidth', 2)
xlabel('Simulations')
ylabel('Price Call')
axis square

hold on
plot(C3, 'linewidth', 2)
legend('MC','Antithetic varaites')
%subplot (2,2,2)
%plot(P1, 'linewidth', 2)
%xlabel('Time')
%ylabel('Price Put')
%axis square

%hold on
%plot(P3, 'linewidth', 2)

subplot (1,2,2)
plot(ErrC1, 'linewidth', 2)
xlabel('Simulations')
ylabel('error Call')
axis square

hold on
plot(ErrC3, 'linewidth', 2)

%subplot (2,2,4)
%plot(ErrP1, 'linewidth', 2)
%xlabel('Time')
%ylabel('error Put')
%axis square

%hold on
%plot(ErrP3, 'linewidth', 2)

legend('MC','Antithetic varaites')