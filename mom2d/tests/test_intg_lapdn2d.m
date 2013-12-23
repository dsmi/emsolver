function test_intg_lapdn2d
% 

rsrc = [ cat(3, [ 1 1 ], [ 2 1 ]); cat(3, [ 1 1 ], [ 1 2 ]); cat(3, [ 0 0 ], [ 1 2 ]) ];
robs = [ cat(3, [ 1 2 ], [ 2 2 ]); cat(3, [ 2 1 ], [ 2 2 ]); cat(3, [ 2 4 ], [ 4 8 ]) ];

vq = intg_lapdn2dq(rsrc,robs);
v = intg_lapdn2d(rsrc,robs);

assertEquals(vq, v, 1e-15)
