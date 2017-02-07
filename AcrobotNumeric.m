clear all
clc

%Robot parameters
m1 = 1/5;
m2 = 1/5;
I1 = 1/20;
I2 = 1/20;
lc1 = 1/10;
lc2 = 1/10;
l1 = 1/5;
l2 = 1/5;


g0 = -981/100;



oldState = [-pi/2 0 0 0]';
newState = [0 0 0 0]';

%0 refers to the desired acceleration
goal = [pi/2 0 0];

Kd = 3 ;
Kp = 2;
deltaT = 15/100;
totalT = 15;
saturationQ2D = 1;
tauLimit = 2;
jointLimitQ1 = pi;
jointLimitQ2 = pi;




totalIterations = ceil(totalT/deltaT);
stateStorage = zeros(totalIterations,2);
stateStorage2 = zeros(totalIterations,6);
taskStorage = zeros(totalIterations,1);
indexStorage = 0;

for t=0:deltaT:totalT
    if (mod(t,1)==0)
        disp(t);
    end
    
     
    tauViolated = false;
    
    indexStorage = indexStorage+1;
    
    
%Robot dynamic parameters
a1 = m1*lc1^2 + I1 + I2 + m2*(l1^2 + lc2^2);
a2 = m2*l1*lc2;
a3 = m2*lc2^2 + I2;

M = [a1 + 2*a2*cos(oldState(2)), a3 + a2*cos(oldState(2));
    a3 + a2*cos(oldState(2)), a3];
C = [-a2*sin(oldState(2))*(oldState(4)^2+2*oldState(3)*oldState(4));
    a2*sin(oldState(2))*oldState(3)^2];

g = g0*[(m1*lc1 + m2*l1)*cos(oldState(1)) + m2*lc2*cos(oldState(1)+oldState(2));
        m2*lc2*cos(oldState(1)+oldState(2))];

    
%Task space definition    
comPos = [(lc1*cos(oldState(1)) + l1*cos(oldState(1)) + lc2*cos(oldState(1)+oldState(2)))/2;
    (lc1*sin(oldState(1)) + l1*sin(oldState(1)) + lc2*sin(oldState(1)+oldState(2)))/2];

comAngle = atan2(comPos(2), comPos(1));


%Jacobian Definition
J = [ 1, (3*cos(oldState(2)) + 1)/(6*cos(oldState(2)) + 10)];
 
Jbar = J(2) - J(1)*inv(M(1,1))*M(1,2);
JbarPinv = pinv(Jbar);

Jdot = [ 0, -(lc2*sin(oldState(2))*(l1 + lc1)*oldState(4)*(l1^2 + 2*l1*lc1 + lc1^2 - lc2^2))/(l1^2 + 2*l1*lc1 + 2*cos(oldState(2))*l1*lc2 + lc1^2 + 2*cos(oldState(2))*lc1*lc2 + lc2^2)^2];

comAngleDot = J*[oldState(3);oldState(4)];

    
   
    
    stateStorage(indexStorage,:) = oldState(1:2);
    taskStorage(indexStorage) = comAngle;
    
    
    
    actualMat = [cos(comAngle), - sin(comAngle); sin(comAngle), cos(comAngle)];
    referenceMat = [cos(goal(1)), - sin(goal(1)); sin(goal(1)), cos(goal(1))];
    errorMat = (actualMat)\referenceMat;
    errorVec = atan2(errorMat(2,1),errorMat(1,1));
    
    vA = goal(3) + Kd*(goal(2) - comAngleDot) + Kp*errorVec;
    
    
    q2DDActual = JbarPinv*(vA -Jdot*[oldState(3);oldState(4)] + J(1)*inv(M(1,1))*(C(1) + g(1)))
   
    
    q1DDActual = -inv(M(1,1))*(M(1,2)*q2DDActual + C(1) + g(1));

    stateStorage2(indexStorage,:) = [oldState', q1DDActual,q2DDActual];
    
    
    tauActual = M(2,1)*q1DDActual + M(2,2)*q2DDActual + C(2) + g(2);

    %inverso = vpa(subs((-M(2,1)*inv(M(1,1))*M(1,2) + M(2,2)),q2,oldState(2)))
    % inverso2 = vpa(subs((-M(2,1)*((M(1,1))\M(1,2)) + M(2,2)),q2,oldState(2)))
    %show2 = vpa(directQ2DDFunc(oldState(1),oldState(2),oldState(3),oldState(4),tauActual))
    
    
     q2DDActual = inv(-M(2,1)*inv(M(1,1))*M(1,2) + M(2,2))*(tauActual + M(2,1)*inv(M(1,1))*(C(1)+g(1))- C(2)-g(2))

     
     q1DDActual = -inv(M(1,1))*(M(1,2)*q2DDActual + C(1) + g(1));
      

    
    newState(3:4) = [oldState(3) + q1DDActual*deltaT;
                    oldState(4) + q2DDActual*deltaT];
                
    if (newState(4)> saturationQ2D)
        newState(4) = saturationQ2D;
        
    elseif (newState(4)< - saturationQ2D)
            newState(4) = - saturationQ2D;
    end
    

    newState(1:2) = [mod(oldState(1) + newState(3)*deltaT + (1/2)*q1DDActual*deltaT^2, 2*pi);
                    mod(oldState(2) + newState(4)*deltaT + (1/2)*q2DDActual*deltaT^2, 2*pi)];

                
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

time_text = text(ax.XLim(1),ax.YLim(1),'');
time_text.FontSize = 14;
time_text.FontWeight = 'bold';

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
    time_text.String = mat2str(i*deltaT);
    
    drawnow;
    
    pause(0.15);
    
    
end


timeStorage = linspace(0,totalT,totalIterations+1)';

figure
title('COM angle evolution');
xlabel('Time');
ylabel('COM angle');
plot(timeStorage,taskStorage);
















