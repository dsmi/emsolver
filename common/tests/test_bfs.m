function test_bfs
%

edges = [ 1 2 ; 2 3 ; 3 4 ; 4 5 ; 5 6 ; 2 5 ];

[ d pred ] = bfs(edges, 1, 6);

test_d = [ 0 ; 1 ; 2 ; 3 ; 2 ; 3 ];
test_pred = [ 1 ; 1 ; 2 ; 3 ; 2 ; 5 ];

assertEquals(test_d, d);
assertEquals(test_pred, pred);
