function testFile2
%testFile2 : run unit test for all file File1
%
%   usage :
%   runTestSuite('testAll');
%   Run several tests, and sum up results. Display the number of test, and
%   the number of passed tests.
%
%   
%   -----
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 31/03/2005
%

%   HISTORY


% first test in this file
a = 3;
assertEquals(a, 3);


newTest;
% second test in this file
a = a+1;
assertEquals(a, 4);

newTest;
% another test. We check approximated value of PI
assertEquals(pi, 3.1415, 1e-3);
assertEquals(pi, 355/113, 1e-3);

newTest
% again another test, with sqrt(2)
b = sqrt(2);
assertEquals(b, 1.4142, 1e-3);


