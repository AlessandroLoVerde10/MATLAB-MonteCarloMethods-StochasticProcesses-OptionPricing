% Linear Congruent Random Number Generator..
clear
clc

a = 5
c = 7
m = 2^8
x0 =0

x_m1=zeros(m,1);

for k=1:m  
    x0=mod((a*x0+c),m);
    x_m1(k) = x0;
end;
disp('The Random Numbers are');
u_m1 = x_m1/m
scatter(u_m1(1:end-1), u_m1(2:end), 'red')
