% MUNIT : unit tests for Matlab
% Version 1.0 - 2005, March 25th.
%
% 
% MUNIT proposes unit tests for matlab.
%
% * Build a test :
%
% Suppose you want to test a function called 'prime', which return prime
% number given by indices. Ex : prime(1) gives 2, prime(3) gives 5,  ...
% You can build a test file called 'testPrime.m', with following code :    
%       assertEquals(prime(1), 2);
%       assertEquals(prime(2), 3);
%       assertEquals(prime(3), 5);
%       ... 
% 
% Two assertions are actually available : 
%       assertTrue(booleanFunction(), boolean);
%       assertEquals(valuedFunction(), expectedValue);
%
% The test is simply run by typing :
%       runTest('testPrime')
%
% If one of the assertion is not satisfied, it will be displayed, as well
% as the file and the line when occuring.
%
%
% * Build a test Suite :
%
% To run several tests together is quite simple :
% write a file 'severalTests.m', with lines :
%       addTest('testFile1');
%       addTest('testFile2');
%       addTestSuite('testOtherFiles');
%       ...
%
% then run it by :
%       runTest('severalTests');
%
% The number and positions of failed assertions is given as well as number
% of passed and ailed tests.
%
%
% -----
%
% File content :
%
% * Create tests :
% assertTrue        assert a variable is true (different to zero)
% assertEquals      assert equality of two variables
% newTest           sum up result current test and start a new one
%
% * Create test suite :
% addTest           add a test in a test suite
% addTestSuite      add a test suite in another test suite
%
% * Run tests :
% runTest           run a matlab unit test
% runTestSuite      run a matlab test suite
%
%
%   -----
%
%   author : David Legland 
%   david.legland@jouy.inra.fr
%   INRA - TPV URPOI - BIA IMASTE
%   created the 13/06/2005.
%

help Contents
