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


%Planning parameters
depthTree = 1;
maxBranching = 1000;
threshold = 0.05;
deltaTPlanning = 0.15;

%Controller parameters
Kd = 1;
Kp = 1;
deltaT = 0.002;
totalSeconds = 0;
totalIterations = totalSeconds/deltaT;

%Constraints on torques and joint positions
tauLimit = 2;
jointLimitQ = pi/12;

%Initial state and task goal state
q = [-pi/2; 0.1; 0.1; 0.1; 0.1];
qD = [0; 0; 0; 0; 0];

goal = [pi/2 0;
        l*5/2 0];
    

cd(tempdir)
file = fullfile(cd, 'JFunc.m');
if exist(file,'file')
    disp('matlabFunctions found... Loading from files...');
else
    disp('matlabFunctions not found... Creating files...')
    createMatlabFunctions(m,l,I,lc,active_joints);
end

primitivesScaling = 2;
primitives = [0.1, 0; -0.1, 0; 0, 0.1; 0, -0.1; 0.1, 0.1; 0.1, -0.1; -0.1, 0.1; -0.1,-0.1; 0 , 0];
primitives = primitives(: ,:);
primitives = primitives*primitivesScaling;

graph = simplePlanning([q', qD'],goal, primitives, tauLimit, jointLimitQ, active_joints, depthTree, maxBranching, threshold, deltaTPlanning, deltaT);


if (size(graph.solutionNodes,2) == 0)
    graph.solutionNodes = size(graph.vectorEdges,2) -200 ;
end



for sol=1:size(graph.solutionNodes,2)

    solutionIndices = graph.vectorEdges{graph.solutionNodes(sol)};


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

    for i=1: size(solutionIndices,2)

        stateStorage = graph.verts(solutionIndices(i),:);

        x1 = [l*cos(stateStorage(1)), l*sin(stateStorage(1))];
        x2 = x1 + [l*cos(stateStorage(1) + stateStorage(2)),l*sin(stateStorage(1) + stateStorage(2))];
        x3 = x2 + [l*cos(stateStorage(1) + stateStorage(2) + stateStorage(3)),l*sin(stateStorage(1) + stateStorage(2)+ stateStorage(3))];
        x4 = x3 + [l*cos(stateStorage(1) + stateStorage(2) + stateStorage(3) + stateStorage(4)),l*sin(stateStorage(1) + stateStorage(2)+ stateStorage(3)+ stateStorage(4))];
        x5 = x4 + [l*cos(stateStorage(1) + stateStorage(2) + stateStorage(3) + stateStorage(4)+ stateStorage(5)),l*sin(stateStorage(1) + stateStorage(2)+ stateStorage(3)+ stateStorage(4)+ stateStorage(5))];

        set(link1,'XData',[x0(1),x1(1)], 'YData',[x0(2),x1(2)])
        set(link2,'XData',[x1(1),x2(1)], 'YData',[x1(2),x2(2)])
        set(link3,'XData',[x2(1),x3(1)], 'YData',[x2(2),x3(2)])
        set(link4,'XData',[x3(1),x4(1)], 'YData',[x3(2),x4(2)])
        set(link5,'XData',[x4(1),x5(1)], 'YData',[x4(2),x5(2)])

        %task_text.String = strcat(mat2str(double(vpa(taskStorage(i,1),3)),3),' COM Angle ',mat2str(double(vpa(taskStorage(i,2),3)),3),' COM Lenght');
        time_text.String = mat2str(i*deltaTPlanning);

        drawnow;

        pause(deltaTPlanning);


    end

end



