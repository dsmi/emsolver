function newTest
%newTest : sum up result current test and start a new one
%
%   usage : newTest;
%   This procedure has to be employed in a script, to differentiate unit
%   tests. Each unit test can be located in a different subfunction, but
%   they can be all in the main function.
%
%   Example :
%
%   function testExample
%   % test first protoype for myFunction
%   res1 = myFunction(4);
%   assertTrue(res1==2);
%   % test second protoype for myFunction
%   newTest;
%   res2 = myFunction(1, 4);
%   assertTrue(res2==1);
%
%   The test program should be lauched using the syntax :
%   runTest('testExample');
%   Then result on the two unit tests will be displayed.
%   
%   
%
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 01/04/2005.
%

if evalin('base', 'exist(''MUNIT'', ''var'');');
    evalin('base', 'MUNIT.tested=MUNIT.tested+1;');
    if evalin('base', 'MUNIT.assertFlag;');
        evalin('base', 'MUNIT.passed=MUNIT.passed+1;');
    end
end
