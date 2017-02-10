function createMatlabFunctions(m,l,I,lc,active_joints)


syms q1 q2 q3 q4 q5 q1D q2D q3D q4D q5D real
n_joints = size(active_joints,1);
omega = [q1D; q1D+q2D; q1D+q2D+q3D; q1D+q2D+q3D+q4D; q1D+q2D+q3D+q4D+q5D];

%Compact definitions
s1 = sin(q1);
s12 = sin(q1+q2);
s123 = sin(q1+q2+q3);
s1234 = sin(q1+q2+q3+q4);
s12345 = sin(q1+q2+q3+q4+q5); 

c1 = cos(q1);
c12 = cos(q1+q2);
c123 = cos(q1+q2+q3);
c1234 = cos(q1+q2+q3+q4);
c12345 = cos(q1+q2+q3+q4+q5); 

%robot link COM
pc1 = [lc*c1;lc*s1];
pc2 = [l*c1+lc*c12;l*s1+lc*s12];
pc3 = [l*c1+l*c12+lc*c123;l*s1+l*s12+lc*s123];
pc4 = [l*c1+l*c12+l*c123+lc*c1234;l*s1+l*s12+l*s123+lc*s1234];
pc5 = [l*c1+l*c12+l*c123+l*c1234+lc*c12345;l*s1+l*s12+l*s123+l*s1234+lc*s12345];

Jc1 = jacobian(pc1,[q1;q2;q3;q4;q5]);
Jc2 = jacobian(pc2,[q1;q2;q3;q4;q5]);
Jc3 = jacobian(pc3,[q1;q2;q3;q4;q5]);
Jc4 = jacobian(pc4,[q1;q2;q3;q4;q5]);
Jc5 = jacobian(pc5,[q1;q2;q3;q4;q5]);

vc1 = Jc1*[q1D;q2D;q3D;q4D;q5D];
vc2 = Jc2*[q1D;q2D;q3D;q4D;q5D];
vc3 = Jc3*[q1D;q2D;q3D;q4D;q5D];
vc4 = Jc4*[q1D;q2D;q3D;q4D;q5D];
vc5 = Jc5*[q1D;q2D;q3D;q4D;q5D];


g = [0;9.81];
U = m*g'*(pc1 + pc2 + pc3 + pc4 + pc5);  


Tv = (1/2)*m*(vc1'*vc1 + vc2'*vc2 + vc3'*vc3 + vc4'*vc4 + vc5'*vc5);
Tomega = simplify((1/2)*I*(omega'*omega));

T = Tv + Tomega;
BSymbols = Bmatrix(T,n_joints);
CSymbols = Cmatrix(BSymbols);
hSymbols = jacobian(U,[q1;q2;q3;q4;q5])';

qSymbols_ordered = reorderingMatrix([q1;q2;q3;q4;q5],active_joints);
qDSymbols_ordered = reorderingMatrix([q1D;q2D;q3D;q4D;q5D],active_joints);

BSymbols = reorderingMatrix(BSymbols, active_joints);
CSymbols = reorderingMatrix(CSymbols, active_joints);
hSymbols = reorderingMatrix(hSymbols, active_joints);

matlabFunction(BSymbols,'File', 'BFunc.m');%@(q2,q3,q4,q5)
matlabFunction(CSymbols,'File', 'CFunc.m');%@(q2,q3,q4,q5,q1D,q2D,q3D,q4D,q5D)
matlabFunction(hSymbols,'File', 'hFunc.m');% @(q1,q2,q3,q4,q5)

comPos = (pc1 + pc2 + pc3 + pc4 + pc5)/n_joints;
comAngle = atan2(comPos(2),comPos(1));
comLength = norm(comPos);
taskSymbols = [comAngle; comLength];
matlabFunction(taskSymbols,'File', 'taskFunc.m');% @(q1,q2,q3,q4,q5)

JSymbols = jacobian(taskSymbols, qSymbols_ordered);
matlabFunction(JSymbols,'File', 'JFunc.m');% @(q1,q2,q3,q4,q5)

Jdq1 = jacobian(JSymbols(:,1),qSymbols_ordered);
Jdq2 = jacobian(JSymbols(:,2),qSymbols_ordered);
Jdq3 = jacobian(JSymbols(:,3),qSymbols_ordered);
Jdq4 = jacobian(JSymbols(:,4),qSymbols_ordered);
Jdq5 = jacobian(JSymbols(:,5),qSymbols_ordered);

JdotSymbols = [Jdq1*qDSymbols_ordered,Jdq2*qDSymbols_ordered,Jdq3*qDSymbols_ordered,Jdq4*qDSymbols_ordered,Jdq5*qDSymbols_ordered];
matlabFunction(JdotSymbols,'File', 'JdotFunc.m');%@(q1,q2,q3,q4,q5,q1D,q2D,q3D,q4D,q5D)

end