 % calcolo del'europea asiatica
clc 
clear
close all
% Model Parameters

S0 = 100;
K = 100;
vol = 0.2;
r = 0.01;

% Nummber of paths Simulated
m1  = 10^1;
m2 = 10^3;
m3  = 10^5;

T = 1;
t = 360;

dt = 1/t;

nsim = 10^2;


% m1 osservazioni

S1 = zeros(m1,t);

S1(:,1) = S0 ;  % Exact Solution

S3 = zeros(m1,t);

S3(:,1) = S0 ; 

Q1 = zeros(m1,t);
Q1(:,1) = S0;
for j = 1: nsim
for i =1:t-1
   dW = sqrt(dt) * randn(1,m1);
    S1(:,i+ 1) = S1(:,i) .* exp( (r-0.5*vol^2)*(dt) + vol*dW');
      S3(:,i+ 1) = S3(:,i) .* exp( (r-0.5*vol^2)*(dt) - vol*dW');
    S = [S3;S1];
end
    Q_T1(:,j) = sqrt(mean(S(:,1:end).^2,2));
    payoffP1 = K-Q_T1;
    payoffC1 = Q_T1 - K;
    P1 = exp(-r*T)*mean(max(payoffP1,0));
    C1 = exp(-r*T)*mean(max(payoffC1,0));
end

PMean1 = mean(P1)
Cmean1 = mean(C1)
Pstdev1 = std(P1)
Cstdev1 = std(C1)
ErrP1 = P1 - PMean1
ErrC1 = C1 - Cmean1

figure
subplot (2,2,1)
plot(C1, 'linewidth', 2)
xlabel('Time')
ylabel('Price Call')
axis square

subplot (2,2,2)
plot(P1, 'linewidth', 2)
xlabel('Time')
ylabel('Price Put')
axis square

subplot (2,2,3)
plot(ErrC1, 'linewidth', 2)
xlabel('Time')
ylabel('Price Call')
axis square

subplot (2,2,4)
plot(ErrP1, 'linewidth', 2)
xlabel('Time')
ylabel('Price Put')
axis square

% 10^3 osservazioni

S2 = zeros(m2,t);

S2(:,1) = S0 ;  %sol. esatta

 S3 = zeros(m2,t);

S3(:,1) = S0 ; 

Q2 = zeros(m2,t);
Q2 (:,1) = S0;
for j = 1: nsim
for i =1:t-1
   dW = sqrt(dt) * randn(1,m2);
    S2(:,i+ 1) = S2(:,i) .* exp( (r-0.5*vol^2)*(dt) + vol*dW');
     S3(:,i+ 1) = S3(:,i) .* exp( (r-0.5*vol^2)*(dt) - vol*dW');
    S = [S3;S2];
   
end
    Q_T2(:,j) = sqrt(mean(S(:,1:end).^2,2));
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

%figure
%subplot (2,2,1)
%plot(C2, 'linewidth', 2)
%xlabel('Time')
%ylabel('Price Call')
%axis square

%subplot (2,2,2)
%plot(P2, 'linewidth', 2)
%xlabel('Time')
%ylabel('Price Put')
%axis square


%subplot (2,2,3)
%plot(ErrC2, 'linewidth', 2)
%xlabel('Time')
%ylabel('Price Call')
%axis square

%subplot (2,2,4)
%plot(ErrP2, 'linewidth', 2)
%xlabel('Time')
%ylabel('Price Put')
%axis square

% 10^5 osservazioni

S3 = zeros(m3,t);

S3(:,1) = S0 ;  %sol. esatta

S2 = zeros(m3,t);

S2(:,1) = S0 ; 

Q3 = zeros(m3,t);
Q3 (:,1) = S0;
for j = 1: nsim
for i =1:t-1
   dW = sqrt(dt) * randn(1,m3);
    S3(:,i+ 1) = S3(:,i) .* exp( (r-0.5*vol^2)*(dt) + vol*dW');
    S2(:,i+ 1) = S2(:,i) .* exp( (r-0.5*vol^2)*(dt) - vol*dW');
    S = [S3;S2];
end
    Q_T3(:,j) = sqrt(mean(S(:,1:end).^2,2));
    payoffP3 = K-Q_T3;
    payoffC3 = Q_T3 - K;
    P3 = exp(-r*T)*mean(max(payoffP3,0));
    C3 = exp(-r*T)*mean(max(payoffC3,0));
end

Pmean3 = mean(P3)
Cmean3 = mean(C3)
Pstdev3 = std(P3)
Cstdev3 = std(C3)

ErrP3 = P3 - Pmean3;
ErrC3 = C3 - Cmean3;

figure
subplot (1,2,1)
plot(C3, 'linewidth', 2)
xlabel('Time')
ylabel('Price Call')
axis square

subplot (1,2,2)
plot(P3, 'linewidth', 2)
xlabel('Time')
ylabel('Price Put')
axis square

subplot (2,2,3)
plot(ErrC3, 'linewidth', 2)
xlabel('Time')
ylabel('Error Call')
axis square

subplot (2,2,4)
plot(ErrP3, 'linewidth', 2)
xlabel('Time')
ylabel('Error Put')
axis square

% plot definitivo

figure
subplot (2,2,1)
plot(C1, 'linewidth', 2)
hold on
plot(C2, 'linewidth', 2)
hold on
plot(C3, 'linewidth', 2)
xlabel('Time')
ylabel('Price Call')
axis square

subplot (2,2,2)
plot(P1, 'linewidth', 2)
hold on
plot(P2, 'linewidth', 2)
hold on
plot(P3, 'linewidth', 2)
xlabel('Time')
ylabel('Price Put')
axis square

subplot (2,2,3)
plot(ErrC1, 'linewidth', 2)
hold on
plot(ErrC2, 'linewidth', 2)
hold on
plot(ErrC3, 'linewidth', 2)
xlabel('Time')
ylabel('Error Call')
axis square

subplot (2,2,4)
plot(ErrP1, 'linewidth', 2)
hold on
plot(ErrP2, 'linewidth', 2)
hold on
plot(ErrP3, 'linewidth', 2)
xlabel('Time')
ylabel('Error Put')
axis square




