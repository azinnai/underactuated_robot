function indicesList = findSortedNeighbors(verts, qRand, jointRanges)

%goalNode = [qRand; 0; 0; 0; 0; 0];
%goalJoints = goalNode(1:5);

%weightsJoints = eye(size(goalJoints,1));  
%weightsJoints = diag([5;5;5;4;3]);
weightTask1 = 0.7;

%goalPotentialEnergy = computePotentialEnergy(goalNode);
goalPotentialEnergy = qRand(2)*sin(qRand(1));


H = zeros(size(verts,1),1);


for i = 1: size(verts,1)

    %nodePotentialEnergy = computePotentialEnergy(verts(i,:));   
    %energyDifference = abs(goalPotentialEnergy - nodePotentialEnergy);
    %H(i) = energyDifference;
    %differenceFromStraightPose = boxMinus(verts(i,1:5)', goalJoints);
    %H(i) = H(i) + norm(differenceFromStraightPose' * weightsJoints * differenceFromStraightPose);
    %H(i) = abs(goalPotentialEnergy - nodePotentialEnergy);

    taskValues = taskFunc(verts(i,1), verts(i,2), verts(i,3),verts(i,4), verts(i,5));
    angleDistance = weightTask1*(boxMinus(qRand(1),taskValues(1)));
    lengthDistance = (1-weightTask1)*(qRand(2) - taskValues(2));
    H(i) = norm([angleDistance; lengthDistance]);

    
end
[Hordered, indicesList] = sort(H,'ascend');



end