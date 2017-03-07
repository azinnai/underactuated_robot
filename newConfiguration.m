function  [qNew , indexArrayUsedPrimitive, state] = newConfiguration( xRand, qNear, xGoal, motionPrimitives, indexArrayFreePrimitives, deltaTPlanning, deltaT, tauLimit, jointLimit, threshold, Knull, active_joints)


nPrimitives = size(motionPrimitives,1);
dist = inf;
simulation = false;

% initialization
qNew = [];
indexArrayUsedPrimitive = [];
indexSuccessfulPrimitive = [];

% loop
for i=1:nPrimitives
    
    if indexArrayFreePrimitives(i)==0 % the i-th motion primive has been already used from qNear
        continue;
    end
    
    % generate the i-th motion primitive
    desiredTaskAcceleration = motionPrimitives(i);
    
    Q = checkConstraints(qNear, desiredTaskAcceleration, Knull, deltaTPlanning, deltaT, tauLimit, jointLimit, active_joints, simulation)

    % now Q represents the motion primitive path -> it is a sequence of configurations
    
    % check path for collisions
    if (Q == 9999) % if the last configuration of the path is not collision free
            indexArrayUsedPrimitive = [ i indexArrayUsedPrimitive]; % the i-th motion primitive brings to collision
            continue; % jump to next cycle 
    end
    
    % compute the distance between the last configuration of the path, i.e. Q(end,:), 
    % and qRand and select the path which brings the robot closer to qRand
    
    taskNewValues = taskFunc(Q(1),Q(2),Q(3),Q(4),Q(5));
    angleDiff = boxMinus(taskNewValues(1),xRand(1));
    lengthDiff = taskNewValues(2) - xRand(2);
    d = norm([angleDiff;lengthDiff]);
    
    if d < dist
        qNew  = Q;
        dist = d;
        indexSuccessfulPrimitive = i; % the i-th motion primitive has been successfully used
    end
    
    
end


indexArrayUsedPrimitive = [ indexSuccessfulPrimitive indexArrayUsedPrimitive];

if  ~isempty(indexSuccessfulPrimitive)
    state = 'advanced';
    nodeTaskValues = taskFunc(qNew(1),qNew(2),qNew(3),qNew(4),qNew(5));
    angleDiff = abs(boxMinus(nodeTaskValues(1), xGoal(1)));
    lengthDiff = abs(nodeTaskValues(2) - xGoal(2));

    if ((angleDiff < threshold(1))&&(lengthDiff < threshold(2)))
        disp('WIN')
        disp(nodeTaskValues')  
        state = 'reached';
    end
else
    state = 'trapped';
end


 
