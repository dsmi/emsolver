function test_rmduptris
%

tri=[1 2 3; 4 5 6; 3 1 2];
t = rmduptris(tri);
assertEquals([ 4 5 6 ], t);
