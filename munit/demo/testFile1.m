function testFile1
%testFile1 : run unit test for all file File1
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
assertTrue(true);
assertEquals(1, 1);


newTest;
% second test in this file
%assertTrue(false);
assertEquals(1, 1);
