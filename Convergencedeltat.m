 % calcolo del'europea asiatica
clc 
clear
close all
%Parametri del modello
vol = 0.1
r = 0.03
S0 = 100
T = 1
m  = 10^4
t = 100
dt = 1/t

nsim = 10^4
K = 98

tt = linspace(0, T, t);
S = zeros(m,t);

S(:,1) = S0 ;  %sol. esatta

Q = zeros(m,t);
Q(:,1) = S0;
for j = 1: nsim
for i =1:t-1
   dW = sqrt(dt) * randn(1,m);
    S(:,i+ 1) = S(:,i) .* exp( (r-0.5*vol^2)*(dt) + vol*dW');
end
    Q_T(:,j) = sqrt(mean(S(:,1:end).^2,2));
    payoffP = K-Q_T;
    payoffC = Q_T - K;
    P = exp(-r*T)*mean(max(payoffP,0));
    C = exp(-r*T)*mean(max(payoffC,0));
end

PMean = mean(P)
Cmean = mean(C)
Pstdev = std(P)
Cstdev = std(C)

figure
subplot (1,2,1)
plot(C, 'linewidth', 2)
xlabel('Time')
ylabel('Price Call')
axis square

subplot (1,2,2)
plot(P, 'linewidth', 2)
xlabel('Time')
ylabel('Price Put')
axis square
