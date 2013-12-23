function test_uncat
%
%

A = rand(3,3,3);
[ A1, A2, A3 ] = uncat(3, A);

assertEquals(A(:,:,1), A1);
assertEquals(A(:,:,2), A2);
assertEquals(A(:,:,3), A3);
