function test_mst
%

edges = [ 1 2 ; 2 3 ; 3 4 ; 4 2 ; 1 3 ];

[ se, sei ] = mst(edges);

sei_test = [ 1 ; 4 ; 5 ];
assertEquals(sei_test, sei);
