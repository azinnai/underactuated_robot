function  simplePlanning(qStart, taskGoal, motionPrimitiveCommandArray, tauLimit, jointLimitQ, active_joints, depthTree, maxBranching, threshold, deltaT)



% Initialization

graph.verts = [];                             % list of vertices: graph.verts(i,:) is the i-th stored configuration  
 
graph.adjMat = sparse(200,200,0); % adjacency matrix of an oriented graph
% graph.adjMat(i,j) = 0 -> no edge exists between between graph.verts(i,:) and graph.verts(j,:)
%                   = d -> 
%                     if (d == 1)
%                        an edge exists from graph.verts(i,:) to graph.verts(j,:)
%                     elseif (d == -1)
%                        an edge exists from graph.verts(j,:) to graph.verts(i,:)
%                     end 
% Hence, if we are building a tree, the parent vertex of graph.verts(i) is the unique j-th vertex such that graph.adjMat(i,j) < 0,
% all the other vertices such that graph.adjMat(i,j) > 0 are children of the i-th vertex. 

graph.adjMat4Edges = sparse(200,200,0); % adjacency matrix for edges identification
graph.vectorEdges = {};                             % list of edges; each edge represents a path in the form of a sequence of configurations [q1; q2; ...; qm]
% graph.adjMat4Edges(i,j) = 0 -> no edge exists between between graph.verts(i,:) and graph.verts(j,:)
%                         = e -> 
%                          if e>0   
%                            graph.vectorEdges(abs(e)) is the oriented path connecting graph.verts(i,:) to graph.verts(j,:) in the C-space manifold 
%                          else
%                            graph.vectorEdges(abs(e)) is the oriented path connecting graph.verts(j,:) to graph.verts(i,:) in the C-space manifold
%                          end
nPrimitives = size(motionPrimitiveCommandArray,1); % number of motion primitives

graph.verts(1,:) = qStart; % insert first vertex


% Main loop 
addedNodes = 0;
search = true;
for s = 1:depthTree
    if (search)
        first = size(graph.verts,1) - addedNodes;
        last = size(graph.verts,1);
        storage = [];
        parents = [];
        for i=first:last
            for j=1:nPrimitives

               newNode = checkConstraints(graph.verts(i,:),motionPrimitiveCommandArray(j,:)', deltaT, tauLimit, jointLimitQ, active_joints);

               if (newNode ~= 9999)
                   storage(end+1,:) = newNode;
                   parents(size(storage,1)) = i;
               end
            end
        end

        for i=1:size(storage,1)
            taskValue = taskFunc(storage(i,1),storage(i,2),storage(i,3),storage(i,4),storage(i,5));
            J = JFunc(storage(i,1),storage(i,2),storage(i,3),storage(i,4),storage(i,5));
            taskDotValue = J* storage(i,6:10)';
            
            dist = norm([taskValue; taskDotValue] - taskGoal');
            if (dist <= threshold)
                disp('WIN!!!!!');%IMPLEMENTARE FUNZIONE DISTANZA
                search = false;
            end
        end

        if (size(storage,1) == 0)
            disp('Failure!!!');
            search = false;
        end
    

        if (size(storage,1) >maxBranching)
             conservedNodes = ones(size(storage,1),1);%DA CANCELLARE
            %PRUNING
        else
            conservedNodes = ones(size(storage,1),1);
        end

        for i=1:size(storage,1)
            if (conservedNodes(i))
                graph.verts(end+1,:) = storage(i,:);
                graph.adjMat(parents(i),size(graph.verts,1)) = 1;
                graph.adjMat(size(graph.verts,1),parents(i)) = -1;
            end
        end

    addedNodes = sum(conservedNodes(:) ==1)
    end

end
    














end