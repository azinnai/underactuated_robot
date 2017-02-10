clear all
clc


%robot parameters for 5r robots with all equal links
m = 0.2;
l = 0.2;
I = 0.05;
lc = 0.1; 

n_joints = 5;
active_joints = [0;1;1;1;1];
n_joints_unactive = sum(active_joints(:) ==0);
n_joints_active = n_joints - n_joints_unactive;
%%%%%%%%%%%%%%%%INITIAL STATE AND GOAL
q = [-pi/2; 0; 0; 0; 0];
qD = [0; 0; 0; 0; 0];

goal = [pi/2 0 0;
        l*5/2 0 0];
    
Kd = 1;
Kp = 1;
deltaT = 0.001;
totalSeconds = 2;
totalIterations = totalSeconds/deltaT;
saturationQD = 5;
tauLimit = 2;
jointLimitQ = pi;

%%%%%%%%%%MATLAB FUNCTION DEFINITION
syms q1 q2 q3 q4 q5 q1D q2D q3D q4D q5D real
omega = [q1D; q1D+q2D; q1D+q2D+q3D; q1D+q2D+q3D+q4D; q1D+q2D+q3D+q4D+q5D];

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

BFunc = matlabFunction(BSymbols);%@(q2,q3,q4,q5)
CFunc = matlabFunction(CSymbols);%@(q2,q3,q4,q5,q1D,q2D,q3D,q4D,q5D)
hFunc = matlabFunction(hSymbols);% @(q1,q2,q3,q4,q5)

comPos = (pc1 + pc2 + pc3 + pc4 + pc5)/n_joints;
comAngle = atan2(comPos(2),comPos(1));
comLength = norm(comPos);
taskSymbols = [comAngle; comLength];
taskFunc = matlabFunction(taskSymbols);% @(q1,q2,q3,q4,q5)

JSymbols = jacobian(taskSymbols, qSymbols_ordered);
JFunc = matlabFunction(JSymbols);% @(q1,q2,q3,q4,q5)

Jdq1 = jacobian(JSymbols(:,1),qSymbols_ordered(1));
Jdq2 = jacobian(JSymbols(:,2),qSymbols_ordered(2));
Jdq3 = jacobian(JSymbols(:,3),qSymbols_ordered(3));
Jdq4 = jacobian(JSymbols(:,4),qSymbols_ordered(4));
Jdq5 = jacobian(JSymbols(:,5),qSymbols_ordered(5));

JdotSymbols = [Jdq1*qDSymbols_ordered(1),Jdq2*qDSymbols_ordered(2),Jdq3*qDSymbols_ordered(3),Jdq4*qDSymbols_ordered(4),Jdq5*qDSymbols_ordered(5) ];
JdotFunc = matlabFunction(JdotSymbols);%@(q1,q2,q3,q4,q5,q1D,q2D,q3D,q4D,q5D)



stateStorage = zeros(totalIterations,n_joints*2);
taskStorage = zeros(totalIterations,4);
errorStorage = zeros(totalIterations,2);
tauStorage = zeros(totalIterations,n_joints_active);

disp('Beginning simulation loop');


%%%%%%%%%%%%%%%%CONTROL LOOP
for t=1:totalIterations
    
    tauViolated = false;
    %Save the state before reordering it
    stateStorage(t,:) = [q',qD'];

    
    B = BFunc(q(2),q(3),q(4),q(5));
    C = CFunc(q(2),q(3),q(4),q(5),qD(1),qD(2),qD(3),qD(4),qD(5));
    h = hFunc(q(1),q(2),q(3),q(4),q(5));
    
    B11 = B(1:n_joints_unactive,1:n_joints_unactive);
    B12 = B(1:n_joints_unactive,n_joints_unactive+1:n_joints);
    B21 = B(n_joints_unactive+1:n_joints,1:n_joints_unactive);
    B22 = B(n_joints_unactive+1:n_joints,n_joints_unactive+1:n_joints);

    C1 = C(1:n_joints_unactive);
    C2 = C(n_joints_unactive+1:n_joints);

    h1 = h(1:n_joints_unactive);
    h2 = h(n_joints_unactive+1:n_joints);


    
    task = taskFunc(q(1),q(2),q(3),q(4),q(5));

    J = JFunc(q(1),q(2),q(3),q(4),q(5));
    
    taskDot = J * qD;

   Jdot = JdotFunc(q(1),q(2),q(3),q(4),q(5),qD(1),qD(2),qD(3),qD(4),qD(5));
     
    J1 = J(:,1:n_joints_unactive);
    J2 = J(:,n_joints_unactive+1:n_joints);    
            
    taskStorage(t,:) = [task' taskDot'];
    
    %Riordino qua, altrimenti avrei problemi ad assegnare parametri alle
    %funzioni
    q = reorderingMatrix(q, active_joints);
    qD = reorderingMatrix(qD, active_joints);

    
    currentAngleMat = [cos(task(1)), - sin(task(1)); sin(task(1)), cos(task(1))];
    referenceAngleMat = [cos(goal(1,1)), - sin(goal(1,1)); sin(goal(1,1)), cos(goal(1,1))];
    errorAngleMat = (currentAngleMat)\referenceAngleMat;
    errorAngleScalar = atan2(errorAngleMat(2,1),errorAngleMat(1,1));
    
    errorVec = [errorAngleScalar;
                (task(2))- goal(2,1)];
            
    errorStorage(t,:) = errorVec';
            
    vA = goal(:,3) + Kd*(goal(:,2) - taskDot) + Kp*errorVec;
 
    %JbarCurrent = J2 -J1*(B11\B12);
    JbarCurrent = J2 - J1*(inv(B11)*B12);
    
    
    JbarPinvCurrent = pinv(JbarCurrent);
    
 
    %q2DDCurrent = JbarPinv * (v - Jdot*qD + J1*(B11\(C1 + h1)));
    q2DDCurrent = JbarPinvCurrent * (vA - Jdot*qD + J1*inv(B11)*(C1 + h1));
    
    %q1DDCurrent = -B11\(B12*q2DD_ordered + C1 + h1);
    q1DDCurrent = -inv(B11)*(B12*q2DDCurrent + C1 + h1);
    
    tauCurrent = B21*q1DDCurrent + B22*q2DDCurrent + C2 + h2;
    
    %show2 = (B22 - B21*(B11\B12))\(tauCurrent + B21*(B11\(C1 + h1)) - C2 - h2);
    %show2 = inv( - B21*inv(B11)*B12 + B22)*(tauCurrent + B21*inv(B11)*(C1 + h1) - C2-h2)
    
    for i=1:size(tauCurrent,1)
        if (tauCurrent(i) > tauLimit)
            tauCurrent(i) = tauLimit;
            tauViolated = true;
        elseif (tauCurrent(i) < - tauLimit)
            tauCurrent(i) = - tauLimit;
            tauViolated = true;
        end
    end

    tauStorage(t,:) = tauCurrent';
    
    if (tauViolated==true)
       disp('tauViolated');
       %q2DDCurrent = (B22 - B21*(B11\B12))\(tauCurrent + B21*(B11\(C1 + h1)) - C2 - h2);
       q2DDCurrent = inv( - B21*inv(B11)*B12 + B22)*(tauCurrent + B21*inv(B11)*(C1 + h1) - C2-h2);
       q1DDCurrent = -inv(B11)*(B12*q2DDCurrent + C1 + h1);

    end


    
    
    
    %%%%%%%%%%%%%%%INTEGRATION
    oldQD = qD;
    qD = qD + [q1DDCurrent; q2DDCurrent]*deltaT;

    for i=1:(n_joints-n_joints_unactive)             
    if (qD(n_joints_unactive + i)> saturationQD)
        qD(n_joints_unactive + i) = saturationQD;
        
    elseif (qD(n_joints_unactive + i)< - saturationQD)
            qD(n_joints_unactive + i) = - saturationQD;
    end
    end
    
    
    q = mod(q + oldQD*deltaT + (1/2)*[q1DDCurrent; q2DDCurrent]*deltaT^2, 2*pi);

     q = INVorder(q,active_joints);
     qD = INVorder(qD,active_joints);

end

  


timeStorage = deltaT:deltaT:totalSeconds;

figure
title('Error Evolution');
xlabel('Time');
ylabel('ERROR');
subplot(2,1,1)
plot(timeStorage,errorStorage(:,1)');
subplot(2,1,2)
plot(timeStorage,errorStorage(:,2)');

figure
title('TAU Evolution');
xlabel('Time');
ylabel('TAU');
legend;
plot(timeStorage,tauStorage);



pause %Wait spacebar press to start real-time simulation 

figure
figureLimits = l*5 + l/2;
xlim([-figureLimits,figureLimits]);
ylim([-figureLimits,figureLimits]);

ax = gca;
task_text = text(ax.XLim(1),ax.YLim(2),'');
task_text.FontSize = 14;
task_text.FontWeight = 'bold';

time_text = text(ax.XLim(1),ax.YLim(1),'');
time_text.FontSize = 14;
time_text.FontWeight = 'bold';

link1 = line('LineWidth',2.5,'Color','r');
if (active_joints(1) == 0)
    link1.Color = 'b';
end
link2 = line('LineWidth',2.5,'Color','r');
if (active_joints(2) == 0)
    link1.Color = 'b';
end
link3 = line('LineWidth',2.5,'Color','r');
if (active_joints(3) == 0)
    link1.Color = 'b';
end
link4 = line('LineWidth',2.5,'Color','r');
if (active_joints(4) == 0)
    link1.Color = 'b';
end
link5 = line('LineWidth',2.5,'Color','r');
if (active_joints(5) == 0)
    link1.Color = 'b';
end


x0 = [0,0]; %Origin of the base link


for i=1:10: size(stateStorage,1)
    x1 = [l*cos(stateStorage(i,1)), l*sin(stateStorage(i,1))];
    x2 = x1 + [l*cos(stateStorage(i,1) + stateStorage(i,2)),l*sin(stateStorage(i,1) + stateStorage(i,2))];
    x3 = x2 + [l*cos(stateStorage(i,1) + stateStorage(i,2) + stateStorage(i,3)),l*sin(stateStorage(i,1) + stateStorage(i,2)+ stateStorage(i,3))];
    x4 = x3 + [l*cos(stateStorage(i,1) + stateStorage(i,2) + stateStorage(i,3) + stateStorage(i,4)),l*sin(stateStorage(i,1) + stateStorage(i,2)+ stateStorage(i,3)+ stateStorage(i,4))];
    x5 = x4 + [l*cos(stateStorage(i,1) + stateStorage(i,2) + stateStorage(i,3) + stateStorage(i,4)+ stateStorage(i,5)),l*sin(stateStorage(i,1) + stateStorage(i,2)+ stateStorage(i,3)+ stateStorage(i,4)+ stateStorage(i,5))];
    
    set(link1,'XData',[x0(1),x1(1)], 'YData',[x0(2),x1(2)])
    set(link2,'XData',[x1(1),x2(1)], 'YData',[x1(2),x2(2)])
    set(link3,'XData',[x2(1),x3(1)], 'YData',[x2(2),x3(2)])
    set(link4,'XData',[x3(1),x4(1)], 'YData',[x3(2),x4(2)])
    set(link5,'XData',[x4(1),x5(1)], 'YData',[x4(2),x5(2)])
    
    task_text.String = strcat(mat2str(double(vpa(taskStorage(i,1),3)),3),' COM Angle ',mat2str(double(vpa(taskStorage(i,2),3)),3),' COM Lenght');
    time_text.String = mat2str(i*deltaT);
    
    drawnow;
    
    %pause(deltaT/2);
    
    
end








