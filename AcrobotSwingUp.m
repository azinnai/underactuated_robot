clear all
clc

syms q1 q2 q1D q2D q1DD q2DD v tau real;

%Robot parameters
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
qDD = [q1DD;q2DD];

g0 = -9.81;

%Robot dynamic parameters
a1 = m1*lc1^2 + I1 + I2 + m2*(l1^2 + lc2^2);
a2 = m2*l1*lc2;
a3 = m2*lc2^2 + I2;

M = [a1 + 2*a2*cos(q2), a3 + a2*cos(q2);
    a3 + a2*cos(q2), a3];
C = [-a2*sin(q2)*(q2D^2+2*q1D*q2D);
    a2*sin(q2)*q1D^2];

g = g0*[(m1*lc1 + m2*l1)*cos(q1) + m2*lc2*cos(q1+q2);
        m2*lc2*cos(q1+q2)];

    
%Task space definition    
comPos = [(lc1*cos(q1) + l1*cos(q1) + lc2*cos(q1+q2))/2;
    (lc1*sin(q1) + l1*sin(q1) + lc2*sin(q1+q2))/2];

comAngle = atan2(comPos(2), comPos(1));


%Jacobian Definition
J = simplify(jacobian(comAngle, q));

Jbar = simplify(J(2) - J(1)*M(1,1)\M(1,2));
JbarPinv = pinv(Jbar);

Jdot = [ 0, -(lc2*sin(q2)*(l1 + lc1)*q2D*(l1^2 + 2*l1*lc1 + lc1^2 - lc2^2))/(l1^2 + 2*l1*lc1 + 2*cos(q2)*l1*lc2 + lc1^2 + 2*cos(q2)*lc1*lc2 + lc2^2)^2];

comAngleDot = J * qD;

%Acceleration and torque equations
requiredQ2DD = JbarPinv * (v -Jdot*qD + J(1)*M(1,1)\(C(1) + g(1)));

requiredQ1DD = -M(1,1)\(M(1,2)*q2DD + C(1) + g(1));

tauCheck = M(2,1)*q1DD + M(2,2)*q2DD + C(2) + g(2);

directQ2DD = (M(2,2)-M(2,1)*M(1,1)\M(1,2))\(tau +M(2,1)*M(1,1)\(C(1)+g(1)) - C(2) - g(2));

state = [q ;qD];
oldState = [-pi/2 0 0 0]';
newState = [0 0 0 0]';

%0 refers to the desired acceleration
goal = [pi/2 0 0];

Kd = 1 ;
Kp = 1;
deltaT = 0.15;
totalT = 0;
saturationQ2D = 1;
tauLimit = 2;
jointLimitQ1 = pi;
jointLimitQ2 = pi;

totalIterations = ceil(totalT/deltaT);
stateStorage = zeros(totalIterations,2);
taskStorage = zeros(totalIterations,1);
indexStorage = 0;

for t=0:deltaT:totalT
    if (mod(t,1)==0)
        disp(t);
    end
    
    tauViolated = false;
    
    indexStorage = indexStorage+1;
    
    taskState = [subs(comAngle,state(1:2), oldState(1:2)) subs(comAngleDot,state, oldState) 0];
  
    
    stateStorage(indexStorage,:) = oldState(1:2);
    taskStorage(indexStorage) = taskState(1);
    
    
    
    
    actualMat = [cos(taskState(1)), - sin(taskState(1)); sin(taskState(1)), cos(taskState(1))];
    referenceMat = [cos(goal(1)), - sin(goal(1)); sin(goal(1)), cos(goal(1))];
    errorMat = (actualMat)\referenceMat;
    errorVec = vpa(atan2(errorMat(2,1),errorMat(1,1)),3);
    
    vA = goal(3) + Kd*(goal(2) - taskState(2)) + Kp*errorVec;

    q2DDActual = subs(requiredQ2DD, [state; v], [oldState; vA]);    
    q1DDActual = subs(requiredQ1DD, [state; q2DD], [oldState; q2DDActual]);
    
    tauActual = vpa(subs(tauCheck,[state;qDD],[oldState;q1DDActual;q2DDActual]),4) 
    
    if (tauActual > tauLimit)
        tauActual
        t
        tauActual = tauLimit;
        tauViolated = true;
    elseif (tauActual < - tauLimit)
        tauActual
        t
        tauActual = - tauLimit;
        tauViolated = true;
    end
    if (tauViolated==true)
      q2DDActual = subs(directQ2DD,[state;tau], [oldState;tauActual]);
      q1DDActual = subs(requiredQ1DD, [state; q2DD], [oldState; q2DDActual]);
     
    end
    
    
    newState(3:4) = [oldState(3) + q1DDActual*deltaT;
                    oldState(4) + q2DDActual*deltaT];
                
    if (newState(4)> saturationQ2D)
        newState(4) = saturationQ2D;
        
    elseif (newState(4)< - saturationQ2D)
            newState(4) = - saturationQ2D;
    end
    
    
  
    
    
    
    
    
    
    
    
    
    
    newState(1:2) = [mod(oldState(1) + newState(3)*deltaT + 0.5*q1DDActual*deltaT^2, 2*pi);
                    mod(oldState(2) + newState(4)*deltaT + 0.5*q2DDActual*deltaT^2, 2*pi)];
    
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
  
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Drawing part

pause %Wait spacebar press to start real-time simulation 

figure
xlim([-0.5,0.5]);
ylim([-0.5,0.5]);

ax = gca;
task_text = text(ax.XLim(1),ax.YLim(2),'');
task_text.FontSize = 14;
task_text.FontWeight = 'bold';

link1 = line;
link1.LineWidth = 2.5;
link1.Color = 'b';
link2 = line;
link2.LineWidth = 2.5;
link2.Color = 'r';

legend([link1,link2], 'UNACTUATED','ACTUATED')

x0 = [0,0]; %Origin of the base link

for i=1:size(stateStorage,1)
    x1 = [l1*cos(stateStorage(i,1)), l2*sin(stateStorage(i,1))];
    x2 = x1 + [l2*cos(stateStorage(i,1) + stateStorage(i,2)),l2*sin(stateStorage(i,1) + stateStorage(i,2))];
    
    set(link1,'XData',[x0(1),x1(1)], 'YData',[x0(2),x1(2)] )
    set(link2,'XData',[x1(1),x2(1)], 'YData',[x1(2),x2(2)])
    
    task_text.String = strcat( mat2str(double(vpa(taskStorage(i),3)),3),' COM Angle');
    
    drawnow;
    
    pause(0.15);
    
    
end
















