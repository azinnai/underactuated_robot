close all
clear all
clc


%robot parameters for 5r robots with all equal links
m = 0.2;
l = 0.2;
I = 0.05;
lc = 0.1; 

active_joints = [0;1;1;1;1];
n_joints = size(active_joints,1);

%Planning parameters
depthTree = 300;
maxBranching = 500;
thresholdAngle = 0.01;
thresholdLength = 0.3;
deltaTPlanning = 0.15;
deltaT = 0.005;

alfa = 0.7;

primitivesScaling = 2;
Knull = 10; %This is used for projected gradient. Should not be a constant. When 0 projected gradient is disabled.

primitives = [1, 0; 
            -1, 0; 
            0, 0.25; 
            0, -0.25; 
            0.5, 0.125; 
            0.5, -0.125; 
            -0.5, 0.125; 
            -0.5,-0.125];
        
primitives = primitives(: ,:);
primitives = primitives*primitivesScaling;

%Constraints on torques and joint positions
tauLimit = 2;
jointLimitQ = pi/12;

%Initial state and task goal state
q = [-pi/2 + 0.1; 0.1; 0.1; 0.1; 0.1];
qD = [0; 0; 0; 0; 0];

goal = [pi/2;l*5/2];

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


threshold = [thresholdAngle; thresholdLength];


graph = simplePlanning([q', qD'],goal, primitives, alfa, Knull, tauLimit, jointLimitQ, active_joints, depthTree, maxBranching, threshold, deltaTPlanning, deltaT);
%graph = RRTplanning([q', qD'],goal, primitives, alfa, Knull, tauLimit, jointLimitQ, active_joints, depthTree, maxBranching, threshold, deltaTPlanning, deltaT);
solutionID = graph.solutionNode;
edgesIDList = graph.vectorEdges{solutionID}';
nodeIterations = ceil(deltaTPlanning/deltaT);

actionsStorage = zeros(size(edgesIDList,1)-1,2);

for i = 1:(size(edgesIDList,1)-1)
    actionsStorage(i,:) = primitives(graph.actionsList{solutionID}(i),:);
end

draw(actionsStorage,active_joints,q,deltaTPlanning,deltaT,Knull);




