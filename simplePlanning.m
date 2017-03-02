function  graph = simplePlanning(qStart, taskGoal, motionPrimitiveArray, Knull, tauLimit, jointLimitQ, active_joints, depthTree, maxBranching, threshold, deltaTPlanning, deltaT)



% Initialization

graph.verts = [];                             % list of vertices: graph.verts(i,:) is the i-th stored configuration  

graph.vectorEdges = {};                             % list of edges; each edge represents a path in the form of a sequence of configurations [q1; q2; ...; qm]

nPrimitives = size(motionPrimitiveArray,1); % number of motion primitives

graph.verts(1,:) = qStart; % insert first vertex
graph.vectorEdges{1} = 1;
graph.actionsList{1} = [];

% Main loop 
search = true;
addedNodes = 0;
numNodesFound = 0;
bestSolution = threshold;
for depth = 1:depthTree

    if (search)
        if (mod(depth,1) == 0)
            disp(strcat('Found nodes: ', mat2str(numNodesFound)));
            disp(strcat('Depth: ', mat2str(depth)));
        end
        first = size(graph.verts,1) - addedNodes;
        if (depth > 1)
            first = first +1;
        end
        last = size(graph.verts,1);
        storage = zeros(maxBranching*nPrimitives,10);
        parents = zeros(maxBranching*nPrimitives,1);
        actionsStorage = zeros(maxBranching*nPrimitives,1);
        numFoundNodes = 0;
        for i=first:last
            for j=1:nPrimitives
               newNode = checkConstraints(graph.verts(i,:),motionPrimitiveArray(j,:)', Knull, deltaTPlanning, deltaT, tauLimit, jointLimitQ, active_joints,false);

               if (newNode ~= 9999)
                   
                   numFoundNodes = numFoundNodes + 1;
                   storage(numFoundNodes,:) = newNode;
                   actionsStorage(numFoundNodes) = j;
                   parents(numFoundNodes) = i;
               end
               
            end

        end
        
        storage = storage(1:numFoundNodes,:);
        parents = parents(1:numFoundNodes);
        actionsStorage = actionsStorage(1:numFoundNodes);
        

        for i=1:numFoundNodes
            
            taskValue = taskFunc(storage(i,1),storage(i,2),storage(i,3),storage(i,4),storage(i,5));
            anglesDiff = boxMinus(taskValue(1), taskGoal(1,1));
            
            if (size(taskGoal,1) == 2)
                lengthDiff = taskValue(2) - taskGoal(2,1);
            else
                lengthDiff = 0;
            end
            
            if (size(taskGoal,2) == 2) %if I have specified velocity goal conditions 
                J = JFunc(storage(i,1),storage(i,2),storage(i,3),storage(i,4),storage(i,5));
                taskDotValue = J* storage(i,6:10)';
                speedDiff = taskDotValue - taskGoal(:,2);
            else speedDiff = 0;
            end
            
            %distanceFromGoal = norm([anglesDiff;lengthDiff;speedDiff]);
            angleDistance = abs(anglesDiff);
            lenghtDistance = abs(lengthDiff);
            
            if (angleDistance <= threshold(1) && lenghtDistance <= threshold(2))
                disp(strcat('WIN!!! at depth ', mat2str(depth)));
                disp(mat2str(storage(i,:)));
                disp(mat2str(taskFunc(storage(i,1),storage(i,2), storage(i,3), storage(i,4), storage(i,5))));
                
                search = false;
                
                newVertexIndex = size(graph.verts,1) + 1;
                
                graph.verts(newVertexIndex,:) = storage(i,:);
                %graph.adjMat(parents(i),newVerticeIndex) = 1;
                %graph.adjMat(newVerticeIndex,parents(i)) = -1;
                
                graph.vectorEdges{newVertexIndex} = [graph.vectorEdges{parents(i)}, newVertexIndex];
                graph.actionsList{newVertexIndex} = [graph.actionsList{parents(i)}, actionsStorage(i)];
                
                if (norm([angleDistance lengthDistance]) <= norm(bestSolution))
                    graph.solutionNode = newVertexIndex;
                    bestSolution = [angleDistance lengthDistance];
                end
                
            end
        end

        if (numFoundNodes == 0)
            disp(strcat('Failure!!! at depth ', mat2str(depth)));
            search = false;
        end
    
        if (search == true) %If I have found a goal, the probabilistic pruning could cut it . So I disable it.
            if (numFoundNodes >maxBranching)
                increment = floor(depth/10);
                alfa = 0.65;
                %alfa = alfa + 0.075 * increment;
                
                conservedNodes = pruningFunction(storage, taskGoal,alfa, maxBranching);
           
            else
                conservedNodes = ones(numFoundNodes,1);
            end

            for i=1:numFoundNodes
                if (conservedNodes(i))
                    newVertexIndex = size(graph.verts,1) + 1;

                    graph.verts(newVertexIndex,:) = storage(i,:);
                    %graph.adjMat(parents(i),newVerticeIndex) = 1;
                    %graph.adjMat(newVerticeIndex,parents(i)) = -1;

                    graph.vectorEdges{newVertexIndex} = [graph.vectorEdges{parents(i)}, newVertexIndex];
                    graph.actionsList{newVertexIndex} = [graph.actionsList{parents(i)}, actionsStorage(i)];
                end
            end
        end

        addedNodes = sum(conservedNodes(:) ==1);
        numNodesFound = size(conservedNodes,1);
    end

end
  


end