function  graph = simplePlanning(qStart, taskGoal, motionPrimitiveArray, tauLimit, jointLimitQ, active_joints, depthTree, maxBranching, threshold, deltaT)



% Initialization

graph.verts = [];                             % list of vertices: graph.verts(i,:) is the i-th stored configuration  

sizeTree = depthTree * maxBranching;

graph.adjMat = sparse(sizeTree,sizeTree,0); % adjacency matrix of an oriented graph
% graph.adjMat(i,j) = 0 -> no edge exists between between graph.verts(i,:) and graph.verts(j,:)
%                   = d -> 
%                     if (d == 1)
%                        an edge exists from graph.verts(i,:) to graph.verts(j,:)
%                     elseif (d == -1)
%                        an edge exists from graph.verts(j,:) to graph.verts(i,:)
%                     end 
% Hence, if we are building a tree, the parent vertex of graph.verts(i) is the unique j-th vertex such that graph.adjMat(i,j) < 0,
% all the other vertices such that graph.adjMat(i,j) > 0 are children of the i-th vertex. 

graph.adjMat4Edges = sparse(sizeTree,sizeTree,0); % adjacency matrix for edges identification
graph.vectorEdges = {};                             % list of edges; each edge represents a path in the form of a sequence of configurations [q1; q2; ...; qm]
% graph.adjMat4Edges(i,j) = 0 -> no edge exists between between graph.verts(i,:) and graph.verts(j,:)
%                         = e -> 
%                          if e>0   
%                            graph.vectorEdges(abs(e)) is the oriented path connecting graph.verts(i,:) to graph.verts(j,:) in the C-space manifold 
%                          else
%                            graph.vectorEdges(abs(e)) is the oriented path connecting graph.verts(j,:) to graph.verts(i,:) in the C-space manifold
%                          end
nPrimitives = size(motionPrimitiveArray,1); % number of motion primitives

graph.verts(1,:) = qStart; % insert first vertex
graph.vectorEdges{1} = 1;

graph.solutionNodes = [];

% Main loop 
addedNodes = 0;
search = true;
for s = 1:depthTree

    if (search)
        s
        first = size(graph.verts,1) - addedNodes;
        last = size(graph.verts,1);
        storage = [];
        parents = [];
        for i=first:last
            for j=1:nPrimitives

               newNode = checkConstraints(graph.verts(i,:),motionPrimitiveArray(j,:)', deltaT, tauLimit, jointLimitQ, active_joints);

               if (newNode ~= 9999)
                   storage(end+1,:) = newNode;
                   parents(end+1) = i;
               end
            end
        end

        for i=1:size(storage,1)
            taskValue = taskFunc(storage(i,1),storage(i,2),storage(i,3),storage(i,4),storage(i,5));
            J = JFunc(storage(i,1),storage(i,2),storage(i,3),storage(i,4),storage(i,5));
            taskDotValue = J* storage(i,6:10)';
            
            angleDiff = boxMinus(taskValue(1), taskGoal(1,1));
            lengthDiff = taskValue(2) - taskGoal(2,1);
            speedDiff = taskDotValue - taskGoal(:,2);
            
            dist = norm([angleDiff;lengthDiff;speedDiff]);
            
            if (dist <= threshold)
                disp(strcat('WIN!!! at depth ', mat2str(s)));
                disp(mat2str(storage(i,:)));
                disp(mat2str(taskFunc(storage(i,1),storage(i,2), storage(i,3), storage(i,4), storage(i,5))));
                
                search = false;
                
                index = size(graph.verts,1) + 1;
                
                graph.verts(index,:) = storage(i,:);
                graph.adjMat(parents(i),index) = 1;
                graph.adjMat(index,parents(i)) = -1;
                
                graph.vectorEdges{index} = [graph.vectorEdges{parents(i)}, index];
                
            end
        end

        if (size(storage,1) == 0)
            disp(strcat('Failure!!! at depth ', mat2str(s)));
            search = false;
        end
    
        if (search == true) %If I have found a goal, the probabilistic pruning could cut it . So I disable it.
            if (size(storage,1) >maxBranching) 
                 conservedNodes = pruningFunction(storage, taskGoal, maxBranching);
            else
                 conservedNodes = ones(size(storage,1),1);
            end

            for i=1:size(storage,1)
                if (conservedNodes(i))
                    index = size(graph.verts,1) + 1;

                    graph.verts(index,:) = storage(i,:);
                    graph.adjMat(parents(i),index) = 1;
                    graph.adjMat(index,parents(i)) = -1;

                    graph.vectorEdges{index} = [graph.vectorEdges{parents(i)}, index];
                end
            end
        end

        addedNodes = sum(conservedNodes(:) ==1);
    end

end
    


end