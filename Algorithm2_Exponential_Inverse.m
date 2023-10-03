%exponential distribution with the inverse method
clear all
close all hidden
clc
lambda = 0.7;
U = rand (10^7,1);
X = -(1/lambda) *  log(U);
Me = mean(X)
Us = sort (U);
Xs = sort (X);
Xd = sort (X, 'descend');
f_x = lambda * (exp(-lambda *X));
f_xd = sort (f_x, 'descend');
F_x = 1 - exp(-lambda*X);
F_xs = sort (F_x);
AX = linspace (0,max(Xs),100000);


figure

 subplot (1,3,1)
plot (Us, Xd, 'linewidth', 2)
title('Exponential random variable');
ylabel ('X')
xlabel ('U')
axis square

subplot (1,3,2)
plot (Xs, f_xd, 'r', 'linewidth', 2)
title('Exponential density');
ylabel ('f(x)')
xlabel ('x')
axis square


subplot (1,3,3)
plot (Xs, F_xs, 'g', 'linewidth', 2)
title ('Exponential CDF')
ylabel ('F(x)')
xlabel ('x')
axis square