function graph = RRTplanning(qStart, xGoal, motionPrimitives, probBiasedExpansion, Knull, tauLimit, jointLimits,  activeJoints,nVertices, dummyParameter, threshold ,deltaTPlanning, deltaT )


nPrimitives = size(motionPrimitives,1); % number of motion primitives

graph.verts = zeros(nVertices,10);                             % list of vertices: graph.verts(i,:) is the i-th stored configuration  
graph.vertsIndexArrayFreePrimitives = zeros(nVertices,nPrimitives);     % graph.vertsIndexArrayFreePrimitives(i,j) = 0 -> the j-th motion primitive has been already used from vertex graph.verts(i,:)    
graph.adjMat = sparse(nVertices,nVertices,0); % adjacency matrix of an oriented graph
graph.vectorEdges = {};                             % list of edges; each edge represents a path in the form of a sequence of configurations [q1; q2; ...; qm]

graph.verts(1,:) = qStart; % insert first vertex
graph.vertsIndexArrayFreePrimitives(1,:) = ones(1,nPrimitives);
graph.vectorEdges{1} = 1;
graph.actionsList{1} = [];
graph.solutionNode = 0;
bestSolutionValue = [inf;inf];



rand('twister', sum(100*clock));  % init random number generator

for s = 1:(nVertices-1)
    
    
    % generate a random configuration
    xRand = randomConfig(jointLimits);    
    
    % goal-biased expansion 
    coin = rand(); % toss a coin
    if coin <= probBiasedExpansion 
        xRand = xGoal; % force xRand = xGoal -> this biases the RRT expansion towards the goal qGoal
    end
    indexqNew = s + 1;

    % Find the configuration in the graph which is closest to xRand and has at least a free motion primitive
    indices = findSortedNeighbors(graph.verts, xRand, s);
    
    for i=1:s
         graph.vertsIndexArrayFreePrimitives(indices(i),:);
        if norm( graph.vertsIndexArrayFreePrimitives(indices(i),:) ) > 0 % if graph.verts(indices(i),:) has at least a free motion primitive 
            qNear = graph.verts(indices(i),:);
            indexqNear = indices(i);
            break;
        end
    end

    
    indexArrayFreePrimitives = graph.vertsIndexArrayFreePrimitives(indexqNear,:); % extract the array representing the free primitives of qNear
    [qNew , indexArrayUsedPrimitive, state] = newConfiguration( xRand, qNear, xGoal, motionPrimitives, indexArrayFreePrimitives, deltaTPlanning, deltaT, tauLimit, jointLimits, threshold, Knull, activeJoints ); 
    %indexArrayUsedPrimitives contains failure primitives and the best available one
    state
    switch state
           
        case 'advanced' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

               taskValues = taskFunc(qNew(1),qNew(2),qNew(3),qNew(4),qNew(5));

               % update the free primitives' array of the vertex qNear
               for h=1:length(indexArrayUsedPrimitive)
                   graph.vertsIndexArrayFreePrimitives(indexqNear,indexArrayUsedPrimitive(h)) = 0;
               end
               
                % insert new vertex in the t.ree structure
               graph.verts(indexqNew,:) = qNew;  
               graph.vertsIndexArrayFreePrimitives(indexqNew,:) = ones(1,nPrimitives);
               
               % insert new edge in the tree structure
               graph.vectorEdges{indexqNew} = [graph.vectorEdges{indexqNear}, indexqNew];
               graph.actionsList{indexqNew} = [graph.actionsList{indexqNear}, indexArrayUsedPrimitive(end)];

               % update adjacency matrices
               graph.adjMat(indexqNear,indexqNew) = 1; 

               if (norm(taskValues - xGoal)<bestSolutionValue)
               	graph.solutionNode = indexqNew;
               	bestSolutionValue = taskValues;
               end
         
        
        case 'reached' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               
               taskValues = taskFunc(qNew(1),qNew(2),qNew(3),qNew(4),qNew(5));

               % update the free primitives' array of the vertex qNear
               for h=1:length(indexArrayUsedPrimitive)
                   graph.vertsIndexArrayFreePrimitives(indexqNear,indexArrayUsedPrimitive(h)) = 0;
               end
               
                % insert new vertex in the tree structure
               graph.verts(indexqNew,:) = qNew;    
               graph.vertsIndexArrayFreePrimitives(indexqNew,:) = ones(1,nPrimitives);
               
               % insert new edge in the tree structure
               graph.vectorEdges{indexqNew} = [graph.vectorEdges{indexqNear}, indexqNew];  
               graph.actionsList{indexqNew} = [graph.actionsList{indexqNear}, indexArrayUsedPrimitive(end)];

               % update adjacency matrices
               graph.adjMat(indexqNear,indexqNew) = 1; 

               graph.solutionNode = indexqNew;
               bestSolutionValue = taskValues;

               disp(indexqNew)

               
               break; % exit from the loop
        
        case 'trapped' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              %All the available primitives failed, sample a new point
        otherwise      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               error('unknown state')
    end    
        

end

graph.adjMat=(graph.adjMat - graph.adjMat'); % the adjacency matrix adjMat is antisymmetric
graph.state = state;


    