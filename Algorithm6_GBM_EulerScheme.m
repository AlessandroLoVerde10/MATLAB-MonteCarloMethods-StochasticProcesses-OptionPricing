
clc 
clear
close all

%Parametri del modello

vol = 0.2;
r = 0.01;
S0 = 100;
T = 1;
m  = 20;
t = 100;
dt = T/t;

tt = linspace(0, T, t);
SE = zeros(m,t);
SE(:,1)= S0 ;  % sol Eulero
Q = zeros(m,t);
Q(:,1) = S0;

for i =1:t-1
   dW = sqrt(dt) * randn(1,m);
   SE(:,i+1) = SE(:,i)+ SE(:,i)*r*(dt) + SE(:,i).*vol.*dW';
end

 
figure

plot(SE', 'linewidth', 2)
xlabel('Time')
ylabel('Price S')
axis square










