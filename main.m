clear all
clc

syms q1 q2 q1dot q2dot real;
syms m1 m2 I1 I2 lc1 lc2 l1 l2 real;

a1 = m1*lc1^2 + I1 + I2 +m2*(l1^2 + lc2^2);
a2 = m2*l1*lc2;
a3 = m2*lc2^2 + I2;

%2R inertia Matrix
M = [a1 + 2*a2*cos(q2), a3 + a2*cos(q2);
    a3 + a2*cos(q2), a3]





