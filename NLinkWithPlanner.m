clear all
clc

%robot parameters for 5r robots with all equal links
m = 0.2;
l = 0.2;
I = 0.05;
lc = 0.1; 


active_joints = [0;1;1;1;1];
n_joints = size(active_joints,1);
n_joints_unactive = sum(active_joints(:) ==0);
n_joints_active = n_joints - n_joints_unactive;


depthTree = 4;
maxBranching = 20;
threshold = 0.25;
deltaTPlanning = 0.15;


q = [-pi/2; 0; 0; 0; 0];
qD = [0; 0; 0; 0; 0];

goal = [pi/2 0;
        l*5/2 0];
    
Kd = 1;
Kp = 1;
deltaT = 0.001;
totalSeconds = 0;
totalIterations = totalSeconds/deltaT;
tauLimit = 2;
jointLimitQ = pi;


createMatlabFunctions(m,l,I,lc,active_joints);

primitives = [0.1, 0; -0.1, 0; 0, 0.1; 0, -0.1; 0.1, 0.1; 0.1, -0.1; -0.1, 0.1; -0.1,-0.1];


simplePlanning([q', qD'],[goal(1,:), goal(2,:)], primitives, tauLimit, jointLimitQ, active_joints,depthTree, maxBranching, threshold, deltaTPlanning);












