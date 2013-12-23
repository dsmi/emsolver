function addTest(testName)
%addTest : add a test in a test suite
%
%   addTest(testName);
%   where testName is a string, will run the corresponding unit test, and
%   add statistics to others tests.
%
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 31/03/2005
%


if evalin('base', 'exist(''MUNIT_SUITE'', ''var'');')
    cmd = sprintf('MUNIT_SUITE.tests = {MUNIT_SUITE.tests{:}, ''%s''};', testName);
    evalin('base', cmd);
end
