function conservedNodes = pruningFunction(newNodes, taskGoal, alfa, maxBranching)


goalNode = [pi/2; 0; 0; 0; 0; 0; 0; 0; 0; 0];
goalJoints = goalNode(1:5);
potentialWeight = 0.5;

%weightsJoints = eye(size(goalJoints,1));  
weightsJoints = diag([5;5;5;4;3]);


[goalPotentialEnergy, goalKineticEnergy] = computeEnergies(goalNode);

%goalPotentialEnergy = taskGoal(2,1)*sin(taskGoal(1,1));
%goalKineticEnergy = 0;



%goalEnergy = potentialWeight * goalPotentialEnergy + (1-potentialWeight) * goalKineticEnergy;

H = zeros(size(newNodes,1),1);
conservedNodes = zeros(size(newNodes,1),1);

for i = 1: size(newNodes,1)
    %taskValue = taskFunc(newNodes(i,1),newNodes(i,2),newNodes(i,3),newNodes(i,4),newNodes(i,5));
    %J = JFunc(newNodes(i,1),newNodes(i,2),newNodes(i,3),newNodes(i,4),newNodes(i,5));
    %taskDotValue = J * newNodes(i,6:10)'; 
    %nodePotentialEnergy = taskValue(2)*sin(taskValue(1));
    %nodeKineticEnergy = 0.5 * (taskDotValue'*taskDotValue);
    
    [nodePotentialEnergy, nodeKineticEnergy] = computeEnergies(newNodes(i,:));
    
    
    %nodeEnergy = potentialWeight*nodePotentialEnergy + (1 - potentialWeight)* nodeKineticEnergy;
    
    energyDifference = potentialWeight * abs(goalPotentialEnergy - nodePotentialEnergy) + (1 - potentialWeight) * abs( goalKineticEnergy - nodeKineticEnergy); 
   
   %energyDifference = abs((nodePotentialEnergy + nodeKineticEnergy) - (goalPotentialEnergy - goalKineticEnergy));
    
    H(i) = 2*(energyDifference);

   % a = H(i)

    
    differenceFromStraightPose = boxMinus(newNodes(i,1:5)', goalJoints);
    
    H(i) = H(i) + (1/7)*norm(differenceFromStraightPose' * weightsJoints * differenceFromStraightPose);

      %  a = H(i)

    
    
end

HBest = min(H);

H = HBest./H;

H = (alfa*H/2) + (1- alfa)*rand(size(H));


Htemp = sort(H,'descend');
Hcut = Htemp(maxBranching + 1);


for i = 1: size(H,1)
    if (H(i)> Hcut)
        conservedNodes(i) = 1;
    else
        conservedNodes(i) = 0;
    end
end



end