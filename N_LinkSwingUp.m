clear all
clc
disp('Computing dynamic model,task terms and compilating matlabFunctions ...')
tic
syms q1 q2 q3 q4 q5 q1D q2D q3D q4D q5D q1DD q2DD q3DD q4DD q5DD tau1 tau2 tau3 tau4 tau5 real;

q = [q1;q2;q3;q4;q5];
qD = [q1D;q2D;q3D;q4D;q5D];
qDD = [q1DD;q2DD;q3DD;q4DD;q5DD];
omega = [q1D; q1D+q2D; q1D+q2D+q3D; q1D+q2D+q3D+q4D; q1D+q2D+q3D+q4D+q5D];

tau = [tau1; tau2; tau3; tau4; tau5];
n_joints = 5;
active_joints = [0;1;0;1;0];
n_joints_unactive = sum(active_joints(:) ==0);
%n_joints_active = n_joints - n_joints_unactive;

%robot parameters for 5r robots with all equal links
m = 1/5;
l = 1/5;
I = 1/20;
lc = 1/10;


%%%%%%%%%%%%%%%%DYNAMIC MODEL COMPUTATION
%Compact definitions
s1 = sin(q1);
s2 = sin(q2);
s3 = sin(q3);
s4 = sin(q4);
s5 = sin(q5);
s12 = sin(q1+q2);
s123 = sin(q1+q2+q3);
s1234 = sin(q1+q2+q3+q4);
s12345 = sin(q1+q2+q3+q4+q5); 

c1 = cos(q1);
c2 = cos(q2);
c3 = cos(q3);
c4 = cos(q4);
c5 = cos(q5);
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

Jc1 = jacobian(pc1,q);
Jc2 = jacobian(pc2,q);
Jc3 = jacobian(pc3,q);
Jc4 = jacobian(pc4,q);
Jc5 = jacobian(pc5,q);

vc1 = Jc1*qD;
vc2 = Jc2*qD;
vc3 = Jc3*qD;
vc4 = Jc4*qD;
vc5 = Jc5*qD;


g = [0;-9.81];
U = m*g'*(pc1 + pc2 + pc3 + pc4 + pc5);  


Tv = (1/2)*m*(vc1'*vc1 + vc2'*vc2 + vc3'*vc3 + vc4'*vc4 + vc5'*vc5);
Tomega = simplify((1/2)*I*(omega'*omega));

T = Tv + Tomega;
B = Bmatrix(T,n_joints);
C = Cmatrix(B);
h = jacobian(U,q)';


B = reorderingMatrix(B, active_joints);
C = reorderingMatrix(C, active_joints);
h = reorderingMatrix(h, active_joints);
tau = reorderingMatrix(tau, active_joints);
tau2 = tau(n_joints_unactive+1:n_joints,1);
q_ordered = reorderingMatrix(q, active_joints);
qD_ordered = reorderingMatrix(qD, active_joints);
qDD_ordered = reorderingMatrix(qDD, active_joints);

q1DD_ordered = qDD_ordered(1:n_joints_unactive);
q2DD_ordered = qDD_ordered(n_joints_unactive+1:n_joints);

q1D_ordered = qD_ordered(1:n_joints_unactive);
q2D_ordered = qD_ordered(n_joints_unactive+1:n_joints);

q1_ordered = q_ordered(1:n_joints_unactive);
q2_ordered = q_ordered(n_joints_unactive+1:n_joints);

B11 = B(1:n_joints_unactive,1:n_joints_unactive);
B12 = B(1:n_joints_unactive,n_joints_unactive+1:n_joints);
B21 = B(n_joints_unactive+1:n_joints,1:n_joints_unactive);
B22 = B(n_joints_unactive+1:n_joints,n_joints_unactive+1:n_joints);

C1 = C(1:n_joints_unactive);
C2 = C(n_joints_unactive+1:n_joints);

h1 = h(1:n_joints_unactive);
h2 = h(n_joints_unactive+1:n_joints);


%%%%%%%%%%%%%%%%TASK TERMS COMPUTATION
comPos = (pc1 + pc2 + pc3 + pc4 + pc5)/n_joints;
comAngle = atan2(comPos(2),comPos(1));
comLength = norm(comPos);
task = [comAngle; comLength]

J = jacobian(task, q_ordered)
%J = simplify(J);

J1 = J(:,1:n_joints_unactive);
J2 = J(:,n_joints_unactive+1:n_joints);


Jdot1 = qD_ordered'*jacobian(J(1,:),q_ordered);
Jdot2 = qD_ordered'*jacobian(J(2,:),q_ordered);
Jdot = [Jdot1;
        Jdot2];

%Jbar = J2 -J1*(B11\B12);
Jbar = J2 - J1*(inv(B11)*B12);
JbarFunc = matlabFunction(Jbar);%%INPUT (q1,q2,q3,q4,q5)


JbarPinv = sym('Jbp',size(Jbar')); %JbarPinv ha dimensioni di Jbar trasposta
v = sym('v',size(task));
taskDot = J * qD;

%requiredQ2DD = JbarPinv * (v - Jdot*qD + J1*(B11\(C1 + h1)));
requiredQ2DD = JbarPinv * (v - Jdot*qD + J1*inv(B11)*(C1 + h1));
requiredQ2DDFunc = matlabFunction(requiredQ2DD); % JbarPinvActual(1,1),JbarPinvActual(1,2),JbarPinvActual(2,1),JbarPinvActual(2,2),oldState(1),oldState(2),oldState(3),oldState(4),oldState(5),oldState(6),oldState(7),oldState(8),oldState(9),oldState(10),vA(1),vA(2))

%requiredQ1DD = -B11\(B12*q2DD_ordered + C1 + h1);
requiredQ1DD = -inv(B11)*(B12*q2DD_ordered + C1 + h1);
requiredQ1DDFunc = matlabFunction(requiredQ1DD);


tauCheck = B21*q1DD_ordered + B22*q2DD_ordered + C2 + h2;
tauCheckFunc = matlabFunction(tauCheck);


%directQ2DD = (B22 - B21*(B11\B12))\(tau2 + B21*(B11\(C1 + h1)) - C2 - h2);
directQ2DD = inv( - B21*inv(B11)*B12 + B22)*(tau2 + B21*inv(B11)*(C1 + h1) - C2-h2);
directQ2DDFunc = matlabFunction(directQ2DD);

matrixCheck = ( - B21*inv(B11)*B12 + B22);
matrixCheckFunc = matlabFunction(matrixCheck);


 
%%%%%%%%%%%%%%%%INITIAL STATE AND GOAL
state = [q_ordered;qD_ordered];
%Se i giunti iniziali sono in [-pi/2 0 0 0 0] c'? una divisione per zero in
%taskDot
oldState = [-pi/2 0.1 0.1 0.1 0 0 0 0 0 0]';
oldState = [reorderingMatrix(oldState(1:n_joints),active_joints); reorderingMatrix(oldState(n_joints+1:n_joints*2),active_joints)];

newState = [0 0 0 0 0 0 0 0 0 0]';

goal = [pi/2 0 0;
        l*5/2 0 0];
    
Kd = 3;
Kp = 1;
deltaT = 15/100;
totalT = 10;
saturationQD = 5;
tauLimit = 2;
jointLimitQ = pi;

totalIterations = ceil(totalT/deltaT);
stateStorage = zeros(totalIterations,n_joints*2);
taskStorage = zeros(totalIterations,size(task,1)*2);
indexStorage = 0;
toc
disp('Beginning simulation loop');

%%%%%%%%%%%%%%%%CONTROL LOOP
for t=0:deltaT:totalT
    tic
    if (mod(t,1)==0)
        disp(t);
    end
    
    tauViolated = false;
    
    indexStorage = indexStorage+1;
    %task state ha dimensione 2x3, come goal. L'ultima colonna inutile serve per fare
    %combaciare le dimensioni.
    taskState = [subs(task,state, oldState) subs(taskDot, state, oldState) [0 0]' ];
  
    
    stateStorage(indexStorage,:) = oldState;
    taskStorage(indexStorage,:) = [taskState(:,1)' taskState(:,2)'];

    
    
    actualAngleMat = [cos(taskState(1,1)), - sin(taskState(1,1)); sin(taskState(1,1)), cos(taskState(1,1))];
    referenceAngleMat = [cos(goal(1,1)), - sin(goal(1,1)); sin(goal(1,1)), cos(goal(1,1))];
    errorAngleMat = (actualAngleMat)\referenceAngleMat;
    errorAngleScalar = atan2(errorAngleMat(2,1),errorAngleMat(1,1));
    
    errorVec = [errorAngleScalar;
                taskState(2,1)- goal(2,1)];
    
    vA = goal(:,3) + Kd*(goal(:,2) - taskState(:,2)) + Kp*errorVec;
 
    JbarActual = JbarFunc(oldState(1),oldState(2),oldState(3),oldState(4),oldState(5));

    JbarPinvActual = pinv(JbarActual);

    q2DDActual = JbarPinvActual * (vA - Jdot*qD + J1*inv(B11)*(C1 + h1));
    %q2DDActual = requiredQ2DDFunc(JbarPinvActual(1,1),JbarPinvActual(1,2),JbarPinvActual(2,1),JbarPinvActual(2,2),oldState(1),oldState(2),oldState(3),oldState(4),oldState(5),oldState(6),oldState(7),oldState(8),oldState(9),oldState(10),vA(1),vA(2));

    q1DDActual = requiredQ1DDFunc(oldState(1),oldState(2),oldState(3),oldState(4),oldState(5),oldState(6),oldState(7),oldState(8),oldState(9),oldState(10),q2DDActual(1),q2DDActual(2));

    tauActual = tauCheckFunc(oldState(1),oldState(2),oldState(3),oldState(4),oldState(5),oldState(6),oldState(7),oldState(8),oldState(9),oldState(10),q1DDActual(1),q1DDActual(2),q1DDActual(3), q2DDActual(1),q2DDActual(2));
    %tauActual = subs(tauCheck,[state;qDD],[oldState;q1DDActual;q2DDActual]); 

    show1 = q2DDActual;
    vpa(show1)
    
    show2 = directQ2DDFunc(oldState(1),oldState(2),oldState(3),oldState(4),oldState(5),oldState(6),oldState(7),oldState(8),oldState(9),oldState(10),tauActual(1),tauActual(2));
    vpa(show2)
    
    show3 = vpa(cond(matrixCheckFunc(oldState(2),oldState(3),oldState(4),oldState(5))))

    for i=1:size(tauActual,1)
        if (tauActual(i) > tauLimit)
            tauActual(i) = tauLimit;
            tauViolated = true;
        elseif (tauActual(i) < - tauLimit)
            tauActual(i) = - tauLimit;
            tauViolated = true;
        end
    end


    if (tauViolated==true)
       disp('tauViolated');
      q2DDActual = directQ2DDFunc(oldState(1),oldState(2),oldState(3),oldState(4),oldState(5),oldState(6),oldState(7),oldState(8),oldState(9),oldState(10),tauActual(1),tauActual(2));  
      q1DDActual = requiredQ1DDFunc(oldState(1),oldState(2),oldState(3),oldState(4),oldState(5),oldState(6),oldState(7),oldState(8),oldState(9),oldState(10),q2DDActual(1),q2DDActual(2));


    end
    
    
   % newState(n_joints+1:n_joints*2) = [oldState(n_joints+1:n_joints+1+n_joints_unactive) + q1DDActual'*deltaT
    %oldState(n_joints+1+n_joints_unactive+1 : n_joints*2) + q2DDActual'*deltaT];

    newState(n_joints+1:n_joints*2) = oldState(n_joints+1:n_joints*2) + [q1DDActual; q2DDActual]*deltaT; 
    
    
   for i=1:(n_joints-n_joints_unactive)             
    if (newState(n_joints + n_joints_unactive + i)> saturationQD)
        newState(n_joints + n_joints_unactive + i) = saturationQD;
        
    elseif (newState(n_joints + n_joints_unactive + i)< - saturationQD)
            newState(n_joints + n_joints_unactive + i) = - saturationQD;
    end
   end
    

    
   % newState(1:n_joints) = [mod(oldState(1:n_joints) + newState(n_joints+1:n_joints+1+n_joints_unactive)*deltaT + 0.5*q1DDActual*deltaT^2, 2*pi);
    %                mod(oldState(2) + newState(4)*deltaT + 0.5*q2DDActual*deltaT^2, 2*pi)];
    newState(1:n_joints) = mod(oldState(1:n_joints) + newState(n_joints+1:n_joints*2)*deltaT + (1/2)*[q1DDActual; q2DDActual]*deltaT^2, 2*pi);
    
    
    
  %  if (newState(1)>jointLimitQ1)
  %      newState(1) = jointLimitQ1;
  %  elseif (newState(1)< -jointLimitQ1)
  %      newState(1) = -jointLimitQ1;
  %  end
       
    
  %  if (newState(2)>jointLimitQ2)
  %      newState(2) = jointLimitQ2;
  %  elseif (newState(2)< -jointLimitQ2)
  %      newState(2) = -jointLimitQ2;
  %  end
                
    oldState = newState
  
    toc

end

  



































