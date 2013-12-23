function [ se, sei ] = mst(edges)
% [ se, sei ] = mst(edges)
% 
% Minimum spanning tree of an unweighted graph.
% Inputs:
%   edges    - nedges-by-2 array of edges
% Outputs:
%   se, sei  - spanning tree edges and indices
%              order of the nodes is not preserved
%

n = max(edges(:)); % Max vertex index
A = sparse(edges(:,1), edges(:,2), 1, n, n);
A = A+A'; % The adjacency matrix

vtree = zeros(n, 1); % vertices in the tree

% search queue
sq = zeros(n, 1);
sqh = 1; % head
sqt = 2; % tail

sq(sqh) = edges(1,1); % start vertex
vtree(edges(1,1)) = 1;

se = [];

while sqt>sqh
    v = sq(sqh); sqh = sqh+1; % pop v off the head of the queue
	va = A(v,:);
	va(find(vtree)) = 0; % drop nodes which are in the tree already
	neigbors = find(va); % neighbor nodes not in the tree
    for w=neigbors
        sq(sqt) = w; sqt = sqt+1; % push to the queue
		vtree(w) = 1;
		se = [ se; [ v w ] ]; % the new spanning edge
    end
end

sei = find(ismember(sort(edges, 2), sort(se, 2), 'rows'));
