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
dt = 1/t;

tt = linspace(0, T, t);
K = 100;
 

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

SEP = Pstdev ;
SEC = Cstdev ;
Pstdev1 = std(P1F);
Cstdev1 = std(C1F);  

ErrP1 = sqrt((P1F - PMean1).^2);
ErrC1 = sqrt((C1F - Cmean1).^2);

%% antitetiche

for m1 = 10^2:10^3

S1 = zeros(m1,t);

S1(:,1) = S0 ;  %sol. esatta

S3 = zeros(m1,t);

S3(:,1) = S0 ; 




for i =1:t-1
   dW = sqrt(dt) * randn(1,m1);
    S1(:,i+ 1) = S1(:,i) .* exp( (r-0.5*vol^2)*(dt) + vol*dW');
      S3(:,i+ 1) = S3(:,i) .* exp( (r-0.5*vol^2)*(dt) - vol*dW');
    S = [S3;S1];
    
    Q_T1 = sqrt(mean(S(:,1:end).^2,2));
    payoffP1 = K-Q_T1;
    payoffC1 = Q_T1 - K;
    P1 = exp(-r*T)*mean(max(payoffP1,0));
    C1 = exp(-r*T)*mean(max(payoffC1,0));
end

P1A(1, m1-10^2+1) = P1;
C1A(1, m1-10^2+1) = C1;

PstdevA (1, m1-10^2+1) = std(P1A);
CstdevA (1, m1-10^2+1) = std(C1A);  
end

PMean1A = mean(P1A)
Cmean1A = mean(C1A)

SEPA = PstdevA;
SECA = CstdevA;
Pstdev1A = std(P1A);
Cstdev1A = std(C1A);  


ErrP1A = sqrt((P1A - PMean1A).^2);
ErrC1A= sqrt((C1A - Cmean1A).^2);



figure
subplot (1,2,1)
plot(C1F, 'linewidth', 1)
xlabel('Number of simulations')
ylabel('Price Call')
axis square

hold on
plot(C1A, 'linewidth', 1)

legend('MC','Antithetic Variates')
%subplot (2,2,2)
%plot(P1F, 'linewidth', 1)
%xlabel('Time')
%ylabel('Price Put')
%axis square

%hold on 
%plot(P1A, 'linewidth', 1)

subplot (1,2,2)
plot(ErrC1, 'linewidth', 1)
xlabel('Number of Simulations')
ylabel('Error Call')
axis square

hold on
plot(ErrC1A, 'linewidth', 1)

hold on
plot (SEC, 'linewidth', 2)

hold on
plot (SECA, 'linewidth', 2)

legend('MC error', 'MC error AV', 'St. error', 'St. error AV')
%subplot (2,2,4)
%plot(ErrP1, 'linewidth', 1)
%xlabel('Time')
%ylabel('Error Put')
%axis square

%hold on
%plot(ErrP1A, 'linewidth', 1)

%hold on
%plot (SEP)

%hold on
%plot (SEPA)
