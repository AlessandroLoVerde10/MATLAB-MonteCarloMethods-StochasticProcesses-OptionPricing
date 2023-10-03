 % calcolo del'europea asiatica
clc 
clear
close all
%Parametri del modello
vol = 0.2;
r = 0.01;
S0 = 100;
T = 1;
t=100;
dt = T/t,

tt = linspace(0, T, t);
K = 100
 

for m1 = 10^2:10^3

S1 = zeros(m1,t);

S1(:,1) = S0 ;  %sol. esatta

 
for i =1:t-1
   dW = sqrt(dt) * randn(1,m1);
    S1(:,i+ 1) = S1(:,i) .* exp( (r-0.5*vol^2)*(dt) + vol*dW');
    
 Q_T1 = sqrt(mean(S1(:,1:end).^2,2));
    payoffP1 = K-Q_T1;
    payoffC1 = Q_T1 - K;
    P1 = exp(-r*T)*mean(max(payoffP1,0));
    C1 = exp(-r*T)*mean(max(payoffC1,0));
end
P1F(1, m1-10^2+1) = P1;
C1F(1, m1-10^2+1) = C1;

Pstdev (1, m1-10^2+1) = std(P1F);
Cstdev (1, m1-10^2+1) = std(C1F);  
end

PMean1 = mean(P1F)
Cmean1 = mean(C1F)

SEP = Pstdev 
SEC = Cstdev 
Pstdev1 = std(P1F);
Cstdev1 = std(C1F);  


ErrP1 = sqrt((P1F - PMean1).^2);
ErrC1 = sqrt((C1F - Cmean1).^2);

figure
subplot (2,2,1)
plot(C1F, 'linewidth', 1)
xlabel('Time')
ylabel('Price Call')
axis square

subplot (2,2,2)
plot(P1F, 'linewidth', 1)
xlabel('Time')
ylabel('Price Put')
axis square

subplot (2,2,3)
plot(ErrC1, 'linewidth', 1)
xlabel('Time')
ylabel('Error Call')
axis square

hold on
plot (SEC)

subplot (2,2,4)
plot(ErrP1, 'linewidth', 1)
xlabel('Time')
ylabel('Error Put')
axis square

hold on
plot (SEP)

