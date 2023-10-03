clc 
clear all
% Typical algorithm for acceptance Rejection Method

N = 10^7;
U1 = rand (N,1);
U2 = rand (N,1);
for i = 1 :N
R(i) = - 2*log(U1(i));
V (i)= 2 * pi* U2(i);
Z1(i) =  sqrt(R(i))*cos(V(i));
Z2(i) =  sqrt(R(i))*sin(V(i));
end
figure 
subplot(2,2,[3,4])
hist3([Z1',Z2'], 'nbins',[50,50])
xlabel('Z1')
ylabel ('Z2')
zlabel ('f(Z1, Z2)')
title ('Bivariate density of (Z1,Z2)')
%figure
%[X,Y] = meshgrid(1:0.5:4,1:4);
%V = 1/(2*pi)*exp((1/2)*Z1'*Z2);
%surf(X,Y,V)

subplot(2,2,1)
histfit(Z1)
xlabel('Z1')
ylabel ('f(Z1)')
title ('Univariate density of Z1')

subplot(2,2,2)
histfit(Z2)
xlabel('Z2')
ylabel ('f(Z2)')
title ('Univariate density of Z2')
%the variable Y has the standard normal distribution