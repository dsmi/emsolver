function [ d pred ] = bfs(edges, s, target)
% [ d pred ] = bfs(edges, s, target)
% 
% Breadth-first search.
% Inputs:
%   e      - nedges-by-2 array of edges
%   s      - starting node
%   target - target node
% Outputs:
%   d      - distances to the nodes, -1 if not reachable
%   pred   - predecessors
%

n = max(edges(:)); % Max vertex index
A = sparse(edges(:,1), edges(:,2), 1, n, n);
A = A+A'; % The adjacency matrix

d=-1*ones(n,1);
pred=zeros(n,1);

% search queue
sq = zeros(n, 1);
sqh = 1; % head
sqt = 2; % tail
sq(sqh) = s; % start from s

d(s)=0;
pred(s)=s;
while sqt>sqh
    v = sq(sqh); sqh = sqh+1; % pop v off the head of the queue
	neigbors = find(A(v,:)); % neighbor nodes
    for w=neigbors
        if d(w) < 0
            sq(sqt) = w; sqt = sqt+1; % push to the queue
			pred(w) = v;
            d(w) = d(v)+1; 
            if w == target
				return
			end
        end
    end
end
