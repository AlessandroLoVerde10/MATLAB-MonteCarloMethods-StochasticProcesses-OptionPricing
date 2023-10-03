
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
S = zeros(m,t);

S(:,1) = S0 ;  %sol. esatta

Q = zeros(m,t);
Q(:,1) = S0;

for i =1:t-1
   dW = sqrt(dt) * randn(1,m);
   S(:,i+ 1) = S(:,i) .* exp( (r-0.5*vol^2)*(dt) + vol*dW');
    Q (:, i+1) = (mean(S(:, 1:i+1),2));
end

 
figure
subplot (1,2,1)
plot(S', 'linewidth', 2)
xlabel('Time')
ylabel('Price S')
axis square

subplot (1,2,2)
plot(Q', 'linewidth', 2)
xlabel('Time')
ylabel('Price Q')
axis square








