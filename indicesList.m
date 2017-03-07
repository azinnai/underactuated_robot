function indicesList = findSortedNeighbors(verts, qRand, jointRanges)

goalNode = [qRand; 0; 0; 0; 0; 0];
goalJoints = goalNode(1:5);

weightsJoints = eye(size(goalJoints,1));  
%weightsJoints = diag([5;5;5;4;3]);


goalPotentialEnergy = computePotentialEnergy(goalNode);

%goalEnergy = potentialWeight * goalPotentialEnergy + (1-potentialWeight) * goalKineticEnergy;

H = zeros(size(verts,1),1);


for i = 1: size(verts,1)

    nodePotentialEnergy = computePotentialEnergy(verts(i,:));
        
    energyDifference = abs(goalPotentialEnergy - nodePotentialEnergy);
    
    H(i) = energyDifference;

   % a = H(i)

    
    differenceFromStraightPose = boxMinus(verts(i,1:5)', goalJoints);
    
    H(i) = H(i) + norm(differenceFromStraightPose' * weightsJoints * differenceFromStraightPose);

      %  a = H(i)
    
end
[Hordered, indicesList] = sort(H,'descend');

end