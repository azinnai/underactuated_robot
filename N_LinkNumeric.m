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
q = [-pi/2; 0.1; 0.1; 0.1; 0.1];
qD = [0; 0; 0; 0; 0];

goal = [pi/2 0 0;
        l*5/2 0 0];
    
Kd = 5;
Kp = 5;
Knull = 5;
deltaT = 0.001;
totalSeconds = 15;
totalIterations = totalSeconds/deltaT;
saturationQD = 5;
tauLimit = 2;
jointLimitQ = pi;
epsilon = 0.000000008; %Used for numerical differentiation

fileName = fullfile('activeJoints.txt');
if exist(fileName,'file')
    fileID = fopen(fileName,'r');
    lastConfiguration = fscanf(fileID, '%f');
    if (lastConfiguration == active_joints)
        disp('matlabFunctions found... Loading from files...');
    else
        disp('matlabFunctions not found... Creating files...')
        createMatlabFunctions(m,l,I,lc,active_joints);
    end
    
else
    disp('matlabFunctions not found... Creating files...')
    createMatlabFunctions(m,l,I,lc,active_joints);
end



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
    
    qOriginal = q;
    
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
  
 
    JbarCurrent = J2 - J1*(inv(B11)*B12);    
    try
    	JbarPinvCurrent = pinv(JbarCurrent);
    catch
    	JbarPinvCurrent = zeros(size(JbarCurrent'));
    end
    

    nullSpaceAcceleration = zeros(n_joints_active,1);
    for p=1:n_joints_active
    
        qPlus = q;
        qPlus(p + n_joints_unactive) = qPlus(p + n_joints_unactive) + epsilon;
        qMinus = q;
        qMinus(p + n_joints_unactive) = qMinus(p + n_joints_unactive) - epsilon;
        
        qPlus = INVorder(qPlus, active_joints);
        qMinus = INVorder(qMinus, active_joints);
        
        JbarPlus = JbarFunc(qPlus(1),qPlus(2),qPlus(3),qPlus(4),qPlus(5));
        JbarMinus = JbarFunc(qMinus(1),qMinus(2),qMinus(3),qMinus(4),qMinus(5));


        HPlus = sqrt(det(JbarPlus*JbarPlus'));
        HMinus = sqrt(det(JbarMinus*JbarMinus'));

        gradHp = (HPlus - HMinus)/(2*epsilon);

        nullSpaceAcceleration(p) = gradHp;
    
    end
    
    projectedGradient = Knull*(eye(n_joints_active) - JbarPinvCurrent*JbarCurrent)*nullSpaceAcceleration

    q2DDCurrent = JbarPinvCurrent * (vA - Jdot*qD + J1*inv(B11)*(C1 + h1)) + projectedGradient;
    
    q1DDCurrent = -inv(B11)*(B12*q2DDCurrent + C1 + h1);
    
    tauCurrent = B21*q1DDCurrent + B22*q2DDCurrent + C2 + h2;

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
       %disp('tauViolated');
       q2DDCurrent = inv( - B21*inv(B11)*B12 + B22)*(tauCurrent + B21*inv(B11)*(C1 + h1) - C2-h2);
       q1DDCurrent = -inv(B11)*(B12*q2DDCurrent + C1 + h1);

    end


    
    
    
    %%%%%%%%%%%%%%%INTEGRATION
    oldQD = qD;
    qD = qD + [q1DDCurrent; q2DDCurrent]*deltaT;

       
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


for i=1:20: size(stateStorage,1)
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








