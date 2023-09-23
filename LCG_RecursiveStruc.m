% Linear Congruent Random Number Generator..
clear
clc
a=5
c = 7
m=2^8
O = clock
x0 = floor(O(6)*10)
x_m1=zeros(m,1);

for k=1:m  
    x0=mod((a*x0+c),m);
    x_m1(k) = x0;
end;

disp('The Random Numbers are');
u_m1 = x_m1/m
me = mean(u_m1)

subplot (2,1,1)
plot(u_m1)
title ('Full period')
ylabel ('u_i')
xlabel ('# of u_i')
 
 a=5
c = 8
m=2^8;
O = clock
x0 = floor(O(6)*10)
x_m1=zeros(m,1);

for k=1:m  
    x0=mod((a*x0+c),m);
    x_m1(k) = x0;
end;
disp('The Random Numbers are');
u_m1 = x_m1/m
subplot (2,1,2)
plot(u_m1,':r', 'Linewidth',2)
title ('Recursive structure')
ylabel ('u_i')
xlabel ('# of u_i')

% end 