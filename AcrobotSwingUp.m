clear all
clc


syms q1 q2 q1D q2D q1DD q2DD v tau real;

%Robot parameters
m1 = 1/5;
m2 = 1/5;
I1 = 1/20;
I2 = 1/20;
lc1 = 1/10;
lc2 = 1/10;
l1 = 1/5;
l2 = 1/5;

q = [q1;q2];
qD = [q1D;q2D];
qDD = [q1DD;q2DD];

g = 981/100;

%Robot dynamic parameters

M = [m1*lc1^2+I1+m2*(l1^2+lc2^2+2*l1*lc2*cos(q2))+I2, m2*(lc2^2+l1*lc2*cos(q2))+I2; 
    m2*(lc2^2+l1*lc2*cos(q2))+I2, m2*lc2^2+I2];

C = [m2*l1*lc2*sin(q2)*(-q2D^2-2*q1D*q2D);
     m2*l1*lc2*sin(q2)*q1D^2];

G = [m1*lc1*g*cos(q1)+m2*g*(lc2*cos(q1+q2)+l1*cos(q1));
    m2*g*lc2*cos(q1+q2)];

pc1 = [lc1*cos(q1);
        lc1*sin(q1)];
pc2 = [l1*cos(q1)+lc2*cos(q1+q2);
        l1*sin(q1)+lc2*sin(q1+q2)];
comPos = (pc1+pc2)/2;

comAngle = atan2(comPos(2),comPos(1));
comAngleFunc = matlabFunction(comAngle,'Vars',q);


J = jacobian(comAngle,q);
J1Func = matlabFunction(J(1),'Vars',q);
J2Func = matlabFunction(J(2),'Vars',q);

comAngleD = J*qD;
comAngleDFunc = matlabFunction(comAngleD,'Vars',[q;qD]);

Jbar = J(2)-J(1)*inv(M(1,1))*M(1,2);
JbarFunc = matlabFunction(Jbar,'Vars',q);
JbarPinv = pinv(Jbar);
%JbarPinvFunc = matlabFunction(JbarPinv,'Vars',q);

Jdot = (jacobian(J)*qD)';

requiredQ2DD = JbarPinv*(v-Jdot*qD+J(1)*inv(M(1,1))*(C(1)+G(1)));
requiredQ2DDFunc = matlabFunction(requiredQ2DD,'Vars',[q;qD;v]);
requiredQ1DD = -inv(M(1,1))*(M(1,2)*q2DD+C(1)+G(1));
requiredQ1DDFunc = matlabFunction(requiredQ1DD,'Vars',[q;qD;q2DD]);

tauCheck = M(2,1)*q1DD+M(2,2)*q2DD+C(2)+G(2);
tauCheckFunc = matlabFunction(tauCheck,'Vars',[q;qD;qDD]);

directQ2DD = inv(-M(2,1)*inv(M(1,1))*M(1,2)+M(2,2))*(tau-M(2,1)*inv(M(1,1))*(C(1)+G(1))+C(2)+G(2));
directQ2DDFunc = matlabFunction(directQ2DD,'Vars',[q;qD;tau]);

state = [-pi/3; 0];
stateD = [0; 0];

goal = [pi/2; 0; 0];

deltaT = 5/1000;
totalIt = 3000;

Kp = 0.5;
Kd = 0.1;
tauLimit = 2.0;
saturationQD = 2.0;

stateStorage = zeros(totalIt,4);
taskStorage = zeros(totalIt,2);

for i=1:totalIt
    
    tauCorrect = true;
    
    stateStorage(i,:) = [state' stateD'];
    
    task = comAngleFunc(state(1),state(2));
    taskD = comAngleDFunc(state(1), state(2),stateD(1),stateD(2));
    
    taskStorage(i,:) = [task taskD];
    
    currentAngleMat = [cos(task), - sin(task); sin(task), cos(task)];
    referenceAngleMat = [cos(goal(1)), - sin(goal(1)); sin(goal(1)), cos(goal(1))];
    errorMat = (currentAngleMat)\referenceAngleMat;
    errorAngle = atan2(errorMat(2,1),errorMat(1,1));
    
    vCurrent = goal(3) + Kd*(goal(2) - taskD) + Kp*errorAngle;
    
    q2DDCurrent = requiredQ2DDFunc(state(1),state(2),stateD(1),stateD(2),vCurrent)
    
    q1DDCurrent = requiredQ1DDFunc(state(1),state(2),stateD(1),stateD(2),q2DDCurrent);
    
    tauCurrent = tauCheckFunc(state(1),state(2),stateD(1),stateD(2),q1DDCurrent,q2DDCurrent)
    
    %check = directQ2DDFunc(state(1),state(2),stateD(1),stateD(2),tauCurrent);
    
    
    if (tauCurrent>tauLimit)
        i
        tauCurrent = tauLimit;
        tauCorrect = false;
    elseif (tauCurrent < - tauLimit)
        i
        tauCurrent = -tauLimit;
        tauCorrect = false;
    end
    
    if (tauCorrect == false)
        q2DDCurrent = directQ2DDFunc(state(1),state(2),stateD(1),stateD(2),tauCurrent);
        q1DDCurrent = requiredQ1DDFunc(state(1),state(2),stateD(1),stateD(2),q2DDCurrent);
    end
    
    oldStateD = stateD;

    stateD = stateD + deltaT*[q1DDCurrent;q2DDCurrent];
    
    if (stateD(2) > saturationQD)
        stateD(2) = saturationQD
    elseif (stateD(2)< -saturationQD)
        stateD(2) = - saturationQD
    end
    

    
    state = mod(state + deltaT*oldStateD + (1/2)*deltaT^2*[q1DDCurrent;q2DDCurrent],2*pi);



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
    x1 = [l1*cos(stateStorage(i,1)), l1*sin(stateStorage(i,1))];
    x2 = x1 + [l2*cos(stateStorage(i,1) + stateStorage(i,2)),l2*sin(stateStorage(i,1) + stateStorage(i,2))];
    
    set(link1,'XData',[x0(1),x1(1)], 'YData',[x0(2),x1(2)] )
    set(link2,'XData',[x1(1),x2(1)], 'YData',[x1(2),x2(2)])
    
    task_text.String = strcat( mat2str(double(vpa(taskStorage(i),3)),3),' COM Angle');
    time_text.String = mat2str(i*deltaT);
    
    drawnow;
    
    pause(deltaT);
    
    
end




%figure
%title('COM angle evolution');
%xlabel('Time');
%ylabel('COM angle');
%plot(timeStorage,taskStorage);
















