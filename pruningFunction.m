function conservedNodes = pruningFunction(newNodes, taskGoal, maxBranching)

alfa = 0.75;
goalJoints = [pi/2; 0 ; 0 ; 0 ; 0];

weightsJoints = diag([5;4;3;2;1]);

goalPotentialEnergy = taskGoal(2,1)*sin(taskGoal(1,1));
goalKineticEnergy = 0;

goalEnergy = goalPotentialEnergy + goalKineticEnergy;

H = zeros(size(newNodes,1),1);
index = 1:size(newNodes,1);
conservedNodes = zeros(size(newNodes,1),1);

for i = 1: size(newNodes,1)
    taskValue = taskFunc(newNodes(i,1),newNodes(i,2),newNodes(i,3),newNodes(i,4),newNodes(i,5));
    J = JFunc(newNodes(i,1),newNodes(i,2),newNodes(i,3),newNodes(i,4),newNodes(i,5));
    taskDotValue = J * newNodes(i,6:10)';
    
    nodePotentialEnergy = taskValue(2)*sin(taskValue(1));
    nodeKineticEnergy = 0.5 * (taskDotValue'*taskDotValue);
    
    nodeEnergy = nodePotentialEnergy + nodeKineticEnergy;
    
    H(i) = abs(goalEnergy - nodeEnergy);
    
    differenceFromStraightPose = boxMinus(newNodes(i,1:5)', goalJoints);
    
    H(i) = H(i) + (1/det(weightsJoints))*norm(differenceFromStraightPose' * weightsJoints * differenceFromStraightPose );
    
    
end

HBest = min(H);

H = HBest./H;

H = (alfa*H/2) + (1- alfa)*rand(size(H));



%for i = 1 : size(H,1)
%    for j = 1 : size(H,1)-1
%        if H(j)<H(j+1)
%            H_temp = H(j+1);
%            index_temp = index(j+1);
%            H(j+1) = H(j);
%            H(j) = H_temp;
%            index(j+1) = index(j);
%            index(j) = index_temp;
%        end
%    end
%end


Htemp = sort(H,'descend');
Hcut = Htemp(maxBranching + 1);


for i = 1: size(H,1)
    if H(i)> Hcut
        conservedNodes(i) = 1;
    else
        conservedNodes(i) = 0;
    end
end



end