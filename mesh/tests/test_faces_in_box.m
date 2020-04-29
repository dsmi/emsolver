function test_faces_in_box
%

[ t1, x1, y1, z1 ] = mksphere(1);
x1 = x1 + 2;

[ t2, x2, y2, z2 ] = mksphere(1);
y2 = y2 + 4;

[ t3, x3, y3, z3 ] = mksphere(1);
z3 = z3 - 6;

[ t x y z ] = joinmeshes({ t1 t2 t3 }, { x1 x2 x3 }, { y1 y2 y3 }, { z1 z2 z3});

% Number of triangles in each of the spheres
nts = size( t1, 1 );

f1 = faces_in_box( t, x, y, z, 0.9, -1.1, -1.1, 3.1, 1.1, 1.1 );
assertEquals( 1:nts, sort( f1 ) )

f2 = faces_in_box( t, x, y, z, -1.1, 2.9, -1.1, 1.1, 5.1, 1.1 );
assertEquals( nts+1:2*nts, sort( f2 ) )

f3 = faces_in_box( t, x, y, z, -1.1, -1.1, -7.1, 1.1, 1.1, -4.9 );
assertEquals( 2*nts+1:3*nts, sort( f3 ) )
