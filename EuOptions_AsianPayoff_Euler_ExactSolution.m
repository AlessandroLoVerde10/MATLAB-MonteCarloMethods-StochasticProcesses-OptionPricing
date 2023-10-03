 % calcolo del'europea asiatica
clc 
clear
close all

%Parametri del modello
vol = 0.2;
r = 0.01;
S0 = 100;
T = 1;
m2 = 10^6;
t = 100;
dt = T/t;
nsim = 10^2;
K = 100;

tt = linspace(0, T, t);

% 10^2 osservazioni

S2 = zeros(m2,t);

S2(:,1) = S0 ;  %sol. esatta

for j = 1: nsim
for i =1:t-1
   dW = sqrt(dt) * randn(1,m2);
    S2(:,i+ 1) = S2(:,i) .* exp( (r-0.5*vol^2)*(dt) + vol*dW');
end
 Q_T2(:,j) = S2(:,end);
    payoffP2 = K-Q_T2;
    payoffC2 = Q_T2 - K;
    P2 = exp(-r*T)*mean(max(payoffP2,0));
    C2 = exp(-r*T)*mean(max(payoffC2,0));
end

Pmean2 = mean(P2)
Cmean2 = mean(C2)
Pstdev2 = std(P2)
Cstdev2 = std(C2)
ErrP2 = P2 - Pmean2;
ErrC2 = C2 - Cmean2;

[Call,Put] = blsprice(S0,K,r,T,vol)

figure
subplot (2,2,1)
plot(C2, 'linewidth', 2)
xlabel('Time')
ylabel('Price Call')
axis square

subplot (2,2,2)
plot(P2, 'linewidth', 2)
xlabel('Time')
ylabel('Price Put')
axis square


subplot (2,2,3)
plot(ErrC2, 'linewidth', 2)
xlabel('Time')
ylabel('Price Call')
axis square

subplot (2,2,4)
plot(ErrP2, 'linewidth', 2)
xlabel('Time')
ylabel('Price Put')
axis square



