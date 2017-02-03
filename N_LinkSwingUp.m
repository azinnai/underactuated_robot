clear all
clc

syms q1 q2 q3 q4 q5 q1D q2D q3D q4D q5D q1DD q2DD q3DD q4DD q5DD tau1 tau2 tau3 tau4 tau5 real;

%syms m l I lc real;

tic

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
m = 0.2;
l = 0.2;
I = 0.05;
lc = 0.1;

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

syms g0 real;

g = [0;g0];
U = -m*g'*(pc1 + pc2 + pc3 + pc4 + pc5);  


Tv = 0.5*m*(vc1'*vc1 + vc2'*vc2 + vc3'*vc3 + vc4'*vc4 + vc5'*vc5);
%Tc = 0.5*m*(vc1'*vc1 + vc2'*vc2);
Tomega = simplify(0.5*I*omega'*omega);

T = Tv + Tomega;
B = Bmatrix(T,n_joints);
C = Cmatrix(B);
h = jacobian(U,q)';


B = reorderingMatrix(B, active_joints);
C = reorderingMatrix(C, active_joints);
h = reorderingMatrix(h, active_joints);
tau = reorderingMatrix(tau, active_joints);
tau(1:n_joints_unactive) = zeros(n_joints_unactive,1);
q_ordered = reorderingMatrix(q, active_joints);
qD_ordered = reorderingMatrix(qD, active_joints);
qDD_ordered = reorderingMatrix(qDD, active_joints);

q1DD_ordered = qDD_ordered(1:n_joints_unactive);
q2DD_ordered = qDD_ordered(n_joints_unactive+1:n_joints);

q1D_ordered = qD_ordered(1:n_joints_unactive);
q2D_ordered = qD_ordered(n_joints_unactive+1:n_joints);

q1_ordered = q_ordered(1:n_joints_unactive);
q2_ordered = q_ordered(n_joints_unactive+1:n_joints);


comPos = (pc1 + pc2 + pc3 + pc4 + pc5)/n_joints;
comAngle = atan2(comPos(2),comPos(1));
comLength = norm(comPos);
task = [comAngle; comLength];

J = jacobian(task, q_ordered);
toc
disp('computed jacobian');
tic
%J = simplify(J);

J1 = J(:,1:n_joints_unactive);
J2 = J(:,n_joints_unactive+1:n_joints);

B11 = B(1:n_joints_unactive,1:n_joints_unactive);
B12 = B(1:n_joints_unactive,n_joints_unactive+1:n_joints);
B21 = B(n_joints_unactive+1:n_joints,1:n_joints_unactive);
B22 = B(n_joints_unactive+1:n_joints,n_joints_unactive+1:n_joints);

C1 = C(1:n_joints_unactive);
C2 = C(n_joints_unactive+1:n_joints);

h1 = h(1:n_joints_unactive);
h2 = h(n_joints_unactive+1:n_joints);



%Jdot = qD_ordered'*jacobian([J(1,:);J(2,:)],q_ordered);
Jdot1 = qD_ordered'*jacobian(J(1,:),q_ordered);
Jdot2 = qD_ordered'*jacobian(J(2,:),q_ordered);
Jdot = [Jdot1;
        Jdot2];

Jbar = J2 -J1*(B11\B12);



%JbarPinv = pinv(Jbar);
%pseudoPt1 = inv(Jbar*Jbar');
%JbarPinv = Jbar'*pseudoPt1;
JbarPinv = sym('Jbp',size(Jbar')); %JbarPinv ha dimensioni di Jbar trasposta
v = sym('v',size(task));
taskDot = J * qD;
toc
disp('dopo pinv');
tic
%%requiredQ2DD = JbarPinv * (v - Jdot*qD + J1*(B11\(C1 + h1)));
requiredQ2DD = JbarPinv * (v - Jdot*qD + J1*inv(B11)*(C1 + h1));
toc
disp('required q2dd');
tic
%requiredQ1DD = -B11\(B12*q2DD_ordered + C1 + h1);
requiredQ1DD = -inv(B11)*(B12*q2DD_ordered + C1 + h1);
toc
disp('required q1dd');
tic

tauCheck = B21*q1DD_ordered + B22*q2DD_ordered + C2 + h2;
toc
disp('taucheck');
tic

%directQ2DD = (B22 - B21*(B11\B12))\(tau(n_joints_unactive+1:n_joints) + B21*(B11\(C1 + h1)) - C2 - h2);
directQ2DD = inv(B22 - B21*(inv(B11)*B12))*(tau(n_joints_unactive+1:n_joints) + B21*(inv(B11)*(C1 + h1)) - C2 - h2);
toc
disp('directq2dd');
tic

state = [q_ordered;qD_ordered];
%Se i giunti iniziali sono in [-pi/2 0 0 0 0] c'? una divisione per zero in
%taskDot
oldState = [-pi/2 0.1 0.1 0.1 0 0 0 0 0 0]';
oldState = [reorderingMatrix(oldState(1:n_joints),active_joints); reorderingMatrix(oldState(n_joints+1:n_joints*2),active_joints)];

newState = [0 0 0 0 0 0 0 0 0 0]';

goal = [pi/2 0 0;
        l*5/2 0 0];
    
Kd = 1;
Kp = 1;
deltaT = 0.15;
totalT = 30;
saturationQD = 1;
tauLimit = 2;
jointLimitQ = pi;

totalIterations = ceil(totalT/deltaT);
stateStorage = zeros(totalIterations,n_joints*2);
taskStorage = zeros(totalIterations,size(task,1)*2);
indexStorage = 0;
toc
disp('inizio simulazione');
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
    toc
    disp('error');
    tic
    JbarActual = subs(Jbar,state,oldState);
    toc
    disp('jbar');
    tic
    JbarPinvActual = pinv(JbarActual);
    toc
    disp('jbarpinv');
    tic
    %Le subs devono stare divise, perch? JbarPinv e v sono simboli 1x1
    %mentre JbarPinvActual e vA sono matrici (di dimensioni diverse) quindi
    %non posso creare un vettore colonna che contiene tutte queste cose.
    %subs NON FUNZIONA! Sostituisce a v ogni elemento di vA, espandendo il
    %vettore iniziale (requiredq2dd). POSSIBILE SOLUZIONE: istanziare v e
    %JbarPinv come matrici simboliche e sostituire element-wise.
    %Non funziona neanche cos?
    for i = 1:size(v,1)
        v(i) = vA(i);
    end
    for i = 1:size(JbarPinv,1)
        for j = 1:size(JbarPinv,2)
            JbarPinv(i,j) =JbarPinvActual(i,j);
        end
    end

    q2DDActual = subs(requiredQ2DD, state, oldState);
    %q2DDActual = subs(requiredQ2DD,{v},{vA});
    %q2DDActual = subs(q2DDActual,{JbarPinv},{JbarPinvActual});

    toc
    disp('q2dd');
    tic
    q1DDActual = subs(requiredQ1DD, [state; q2DD], [oldState; q2DDActual]);
    toc
    disp('q1dd');
    tic
    tauActual = subs(tauCheck,[state;qDD],[oldState;q1DDActual;q2DDActual]); 
    toc
    disp('tau');
    tic
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
      q2DDActual = subs(directQ2DD,[state;tau], [oldState;tauActual]);
      q1DDActual = subs(requiredQ1DD, [state; q2DD], [oldState; q2DDActual]);

    end
    
    
   % newState(n_joints+1:n_joints*2) = [oldState(n_joints+1:n_joints+1+n_joints_unactive) + q1DDActual'*deltaT
    %oldState(n_joints+1+n_joints_unactive+1 : n_joints*2) + q2DDActual'*deltaT];
    newState(n_joints+1:n_joints*2) = oldState(n_joints+1:n_joints*2) + [q1DDActual' q2DDActual']*deltaT; 
    
    
                
    %if (newState(4)> saturationQ2D)
    %    newState(4) = saturationQ2D;
    %    
    %elseif (newState(4)< - saturationQ2D)
    %        newState(4) = - saturationQ2D;
    %end
    

    
   % newState(1:n_joints) = [mod(oldState(1:n_joints) + newState(n_joints+1:n_joints+1+n_joints_unactive)*deltaT + 0.5*q1DDActual*deltaT^2, 2*pi);
    %                mod(oldState(2) + newState(4)*deltaT + 0.5*q2DDActual*deltaT^2, 2*pi)];
    newState(1:n_joints) = mod(oldState(1:n_joints) + newState(n_joints+1:n_joints*2)*deltaT + 0.5*[q1DDActual' q2DDActual']*deltaT^2, 2*pi);
    
    
    
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
                
    oldState = newState;
  
    toc
    disp('fine simulazione');
end

  



































