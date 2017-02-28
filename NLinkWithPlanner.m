close all
clear all
clc


%robot parameters for 5r robots with all equal links
m = 0.2;
l = 0.2;
I = 0.05;
lc = 0.1; 

active_joints = [0;1;0;1;1];
n_joints = size(active_joints,1);

%Planning parameters
depthTree = 60;
maxBranching = 3000;
threshold = 0.05;
deltaTPlanning = 0.15;
deltaT = 0.015;
for r = 0
primitivesScaling = 0.3;
Knull = r; %This is used for projected gradient. Should not be a constant. When 0 projected gradient is disabled.

primitives = [1, 0; -1, 0; 0, 1; 0, -1; 1, 1; 1, -1; -1, 1; -1,-1];
primitives = primitives(: ,:);
primitives = primitives*primitivesScaling;
found = false;

%Constraints on torques and joint positions
tauLimit = 2;
jointLimitQ = pi/12;

%Initial state and task goal state
q = [-pi/2 + 0.1; 0.1; 0.1; 0.1; 0.1];
qD = [0; 0; 0; 0; 0];

goal = pi/2;

fileName = 'activeJoints.txt';
if exist(fileName,'file')
    fileID = fopen(fileName,'r');
    lastConfiguration = fscanf(fileID, '%f');
    if (lastConfiguration == active_joints)
        disp('matlabFunctions for the current configuration found... Loading from files...');
    else
        disp('matlabFunctions not found... Creating files...')
        createMatlabFunctions(m,l,I,lc,active_joints);
    end
    
else
    disp('matlabFunctions not found... Creating files...')
    createMatlabFunctions(m,l,I,lc,active_joints);
end





graph = simplePlanning([q', qD'],goal, primitives,Knull, tauLimit, jointLimitQ, active_joints, depthTree, maxBranching, threshold, deltaTPlanning, deltaT);

disp('after planning')
%If I have not found a solution, plot something anyhow.
if sum(strcmp(fieldnames(graph), 'solutionNode')) == 0
   graph.solutionNode = size(graph.vectorEdges,2)  ;
else
   fname = strcat('solution_',mat2str(r), '.mat');
   found = true;
end
solutionID = graph.solutionNode;
edgesIDList = graph.vectorEdges{solutionID}';
nodeIterations = ceil(deltaTPlanning/deltaT);


stateStorage = zeros(size(edgesIDList,1), n_joints*2);
actionsStorage = zeros(size(edgesIDList,1)-1,2);
stateStorageNew = zeros(nodeIterations*size(actionsStorage,1) +1,n_joints*2);
actionsStorageNew = zeros (nodeIterations*size(actionsStorage,1),2);
taskStorage = zeros(size(stateStorageNew,1), 2);
errorStorage = zeros(size(stateStorageNew,1), 1);

for i = 1:(size(edgesIDList,1)-1)
    actionsStorage(i,:) = primitives(graph.actionsList{solutionID}(i),:);
    actionsStorageNew(((i-1)*nodeIterations + 1): i*nodeIterations,:) = repmat(actionsStorage(i,:),nodeIterations,1);
end

for i = 1:(size(edgesIDList,1))
    stateStorage(i,:) = graph.verts(edgesIDList(i),:);    
end



for i = 1:size(stateStorageNew,1)
    if (i == 1)
        stateStorageNew(i,:) =  [q',qD'];
    else
        stateStorageNew(i,:) = checkConstraints(stateStorageNew(i-1,:), actionsStorageNew(i-1,:)', r, deltaT,deltaT,tauLimit,jointLimitQ,active_joints,true);
    end
    
    taskState = taskFunc(stateStorageNew(i,1),stateStorageNew(i,2),stateStorageNew(i,3),stateStorageNew(i,4),stateStorageNew(i,5));

    taskStorage(i,:) = taskState';
   
    currentAngleMat = [cos(taskState(1)), - sin(taskState(1)); sin(taskState(1)), cos(taskState(1))];
    referenceAngleMat = [cos(goal(1,1)), - sin(goal(1,1)); sin(goal(1,1)), cos(goal(1,1))];
    errorAngleMat = (currentAngleMat)\referenceAngleMat;
    errorAngleScalar = atan2(errorAngleMat(2,1),errorAngleMat(1,1));
            
    errorStorage(i,:) = errorAngleScalar;

end




timeStorage = 0:deltaT:(deltaT*(size(stateStorageNew,1) - 1));
errorStorage = abs(errorStorage);
figure
plot(timeStorage,errorStorage(:,1)');
title('Error Evolution')
xlabel('Time')
ylabel('Error')

figure

subplot(2,1,1)
plot(timeStorage,taskStorage(:,1)');
xlabel('Time')
ylabel('Angle')
title('Task 1 Evolution')
subplot(2,1,2)
plot(timeStorage,taskStorage(:,2)');
xlabel('Time')
ylabel('Length')
title('Task 2 Evolution')
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
    link2.Color = 'b';
end
link3 = line('LineWidth',2.5,'Color','r');
if (active_joints(3) == 0)
    link3.Color = 'b';
end
link4 = line('LineWidth',2.5,'Color','r');
if (active_joints(4) == 0)
    link4.Color = 'b';
end
link5 = line('LineWidth',2.5,'Color','r');
if (active_joints(5) == 0)
    link5.Color = 'b';
end


x0 = [0,0]; %Origin of the base link

for i=1: size(stateStorageNew,1)

    x1 = [l*cos(stateStorageNew(i,1)), l*sin(stateStorageNew(i,1))];
    x2 = x1 + [l*cos(stateStorageNew(i,1) + stateStorageNew(i,2)),l*sin(stateStorageNew(i,1) + stateStorageNew(i,2))];
    x3 = x2 + [l*cos(stateStorageNew(i,1) + stateStorageNew(i,2) + stateStorageNew(i,3)),l*sin(stateStorageNew(i,1) + stateStorageNew(i,2)+ stateStorageNew(i,3))];
    x4 = x3 + [l*cos(stateStorageNew(i,1) + stateStorageNew(i,2) + stateStorageNew(i,3) + stateStorageNew(i,4)),l*sin(stateStorageNew(i,1) + stateStorageNew(i,2)+ stateStorageNew(i,3)+ stateStorageNew(i,4))];
    x5 = x4 + [l*cos(stateStorageNew(i,1) + stateStorageNew(i,2) + stateStorageNew(i,3) + stateStorageNew(i,4)+ stateStorageNew(i,5)),l*sin(stateStorageNew(i,1) + stateStorageNew(i,2)+ stateStorageNew(i,3)+ stateStorageNew(i,4)+ stateStorageNew(i,5))];

    set(link1,'XData',[x0(1),x1(1)], 'YData',[x0(2),x1(2)])
    set(link2,'XData',[x1(1),x2(1)], 'YData',[x1(2),x2(2)])
    set(link3,'XData',[x2(1),x3(1)], 'YData',[x2(2),x3(2)])
    set(link4,'XData',[x3(1),x4(1)], 'YData',[x3(2),x4(2)])
    set(link5,'XData',[x4(1),x5(1)], 'YData',[x4(2),x5(2)])

    task_text.String = strcat(mat2str(double(vpa(taskStorage(i,1),3)),3),' COM Angle ',mat2str(double(vpa(taskStorage(i,2),3)),3),' COM Lenght');
    time_text.String = mat2str(i*deltaT);

    drawnow;

end


end

