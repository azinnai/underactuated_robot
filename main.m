clear all
clc

syms q1 q2 q1D q2D q1DD q2DD real;

m1 = 0.2;
m2 = 0.2;
I1 = 0.05;
I2 = 0.05;
lc1 = 0.1;
lc2 = 0.1;
l1 = 0.2;
l2 = 0.2;

q = [q1;q2];
qD = [q1D;q2D];

a1 = m1*lc1^2 + I1 + I2 + m2*(l1^2 + lc2^2);
a2 = m2*l1*lc2;
a3 = m2*lc2^2 + I2;

h = -m2*l1*lc2*sin(q2);

g0 = -9.81;

%2R inertia Matrix symbolic
M = [a1 + 2*a2*cos(q2), a3 + a2*cos(q2);
    a3 + a2*cos(q2), a3];
C = [-a2*sin(q2)*(q2D^2+2*q1D*q2D);
    a2*sin(q2)*q1D^2];

g = g0*[(m1*lc1 + m2*l1)*cos(q1) + m2*lc2*cos(q1+q2);
        m2*lc2*cos(q1+q2)];

comPos = [(lc1*cos(q1) + l1*cos(q1) + lc2*cos(q1+q2))/2;
    (lc1*sin(q1) + l1*sin(q1) + lc2*sin(q1+q2))/2];

comAngle = atan2(comPos(2), comPos(1));

J = simplify(jacobian(comAngle, q));

Jbar = simplify(J(2) - J(1)*M(1,1)\M(1,2));
JbarPinv = pinv(Jbar);

Jdot = [ 0, -(lc2*sin(q2)*(l1 + lc1)*q2D*(l1^2 + 2*l1*lc1 + lc1^2 - lc2^2))/(l1^2 + 2*l1*lc1 + 2*cos(q2)*l1*lc2 + lc1^2 + 2*cos(q2)*lc1*lc2 + lc2^2)^2];

comAngleDot = J * qD;

syms v;

requiredQ2DD = JbarPinv * (v -Jdot*qD + J(1)*M(1,1)\(C(1) + g(1)));

requiredQ1DD = -M(1,1)\(M(1,2)*q2DD + C(1) + g(1));

state = [q ;qD];
oldState = [0 -pi/2 0 0]';
newState = [0 0 0 0]';

%0 refers to the desired acceleration

Kd = 0.1;
Kp = 0.5;
goal = [pi/2 0 0];
deltaT = 0.15;

for t=0:deltaT:10
    
    oldTaskState = [subs(comAngle,state(1:2), oldState(1:2)) subs(comAngleDot,state, oldState) 0];
    
    vpa(oldTaskState,3)
    
    vA = goal(3) + Kd*(goal(2) - oldTaskState(2)) + Kp*(goal(1) - oldTaskState(1));
    
    q2DDActual = subs(requiredQ2DD, [state; v], [oldState; vA]);
    
    q1DDActual = subs(requiredQ1DD, [state; q2DD], [oldState; q2DDActual]);
    
    newState(3:4) = [oldState(3) + q1DDActual*deltaT;
                    oldState(4) + q2DDActual*deltaT];
    
    newState(1:2) = [mod(oldState(1) + newState(3)*deltaT + 0.5*q1DDActual*deltaT^2, 2*pi);
                    mod(oldState(2) + newState(4)*deltaT + 0.5*q2DDActual*deltaT^2, 2*pi)];
    
    oldState = newState
    
end
    














