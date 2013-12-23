function test_shortgndz
%
%

% Test the simplest two-port to one-port case
Z2 = [ -30.3-j*1291.7i -30.3-j*1306.7; -30.3-j*1306.7 -30.3-j*1291.7 ];

Ztest = (Z2(1,1)-Z2(1,2))+(Z2(2,2)-Z2(1,2));
Z = shortgndz(Z2);

assertEquals(Ztest, Z);
