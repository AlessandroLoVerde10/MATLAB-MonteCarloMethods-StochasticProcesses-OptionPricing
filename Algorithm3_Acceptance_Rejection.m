clc 
clear all
% Typical algorithm for acceptance Rejection Method
N = 10^7;
c = 1.3155;
U1 = rand (N,1);
U2 = rand (N,1);
U3 = rand (N,1);
Y = - log(U1);
Ys = sort (Y);
g_y =  c*(exp(-Ys))/2;

Y_0 = zeros(N,1);
for i = 1:N  
Bou(i,1) = exp(-0.5*(Y(i,1)-1)^2);
if U2(i,1) <= Bou (i,1)
    Y_0(i,1) = Y (i,1);
end

if U3(i,1) <= 0.5 %i decide the sign randomly
 Y_2(i,1) = -Y_0(i,1)  ; 
else
  Y_2 (i,1) = Y_0(i,1);
end  
end
X=Y_2;
X_s = sort (X);

f_x =exp(-(X_s.^2)/2)/ sqrt((2)*(3.1415));

figure
plot (X_s, f_x, 'r', 'linewidth', 3)
hold on
plot (-Ys, g_y, 'b', 'linewidth', 3)
hold on
plot (Ys, g_y, 'b', 'linewidth', 3)
xlabel('X, Y')
ylabel('f(X), g(Y)')
legend('f(X) = Normal','g(Y) = Double Exponential')
title ('Normal r.v. X from a Scaled Double Exponential')
axis square
%the variable Y has the distributon of X with density f