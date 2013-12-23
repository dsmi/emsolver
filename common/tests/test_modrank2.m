function test_modrank2
%

A = [ 0 0 0 0 1 ; 0 0 0 1 0 ; 0 0 1 1 0 ];

[ R, x ] = modrank2(A);
assertEquals(3, R);

Ax = [ A; x' ];
Rax = modrank2(Ax);
assertEquals(4, Rax);
