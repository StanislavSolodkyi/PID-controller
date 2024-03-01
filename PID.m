clc;
clear;
close all;
K=2.18
T1 = 0.65
T2 =  1.39
ID = 0.5
G = tf([K/(T1*T2)], [1, (T1+T2)/(T1*T2), 1/(T1*T2)], 'InputDelay', ID)

step(G)
title('povodna prenosova funkcia')


%%optimalny modul
p1 = T1 + T2 + ID
p2 = (T1.^2)+(T2.^2)
p3 = (T1.^3)+(T2.^3)
p4 = (T1.^4)+(T2.^4)
p5 = (T1.^5)+(T2.^5)

a0 = 1
a1 = p1
a2 = ((p1.^2) - p2)/2
a3 = ((p1.^3) - (3*p1*p2) + (2*p3))/6
a4 = ((p1.^4) - (6*(p1.^2)*p2) + (8*p1*p3) + (3 *(p2.^2)) - 6*p4)/24
a5 = ((p1.^5) - (10*(p1.^3)*p2) + (20*(p1.^2)*p3) +(15*p1*(p2.^2)) - (30*p1*p4) - (20*p2*p3) + (24*p5))/120

A = [a1 -1 0; a3 -a2 a1; a5 -a4 a3]
b = (1/(2*K)).*[1; (-(a1.^2)+(2*a2)); ((a2.^2) - (2*a1*a3) + (2*a4) )] 
x = A\b
I_om = x(1)
P_om = x(2)
D_om = x(3)

Gr_om = tf([D_om P_om I_om], [1 0])

Mr1 = feedback(Gr_om * G,1)

figure
hold on
step(Mr1)


%%inverzna dynamika


a = 1/(1.944 * ID)
Ti = T1 + T2
P_inv = (a/K)*Ti

Td = (T1*T2)/(T1+T2)
D_inv = P_inv * Td
I_inv = P_inv/Ti

Gr2 = tf([D_inv P_inv I_inv], [1 0])
Mr2 = feedback(Gr2*G, 1)
step(Mr2)

%metoda casovych konstant
Ts = T1 + T2 + ID
P_sck = 2/K

Ti_sck = 0.8 * Ts
Td_sck = 0.1944 * Ts 

I_sck = P_sck / Ti_sck
D_sck = P_sck * Td_sck

Gr3 = tf([D_sck P_sck I_sck ], [1 0])

Mr3 = feedback(Gr3*G, 1)



step(Mr3)

s1 = sim('Zad1_.slx')
figure
hold on 
plot(s1.s1y)
plot(s1.s1e)
figure 
plot(s1.s1u)

Gd = pade(G,1)
syms x;
[num, den] = tfdata(1/(1+Gd*Gr3))
num = vpa(num)
den = vpa(den)
num = poly2sym(num)
den = poly2sym(den)
e =  limit(num/den, x , 0)

[num, den] = tfdata(Gd*Gr3/(1+Gd*Gr3))
num = vpa(num)
den = vpa(den)
num = poly2sym(num)
den = poly2sym(den)
y = limit(num/den, x , 0)

[num, den] = tfdata(Gr3/(1+Gd*Gr3))
num = vpa(num)
den = vpa(den)
num = poly2sym(num)
den = poly2sym(den)
u = limit(num/den, x , 0)

figure 
nyquist(Gr3*G)



