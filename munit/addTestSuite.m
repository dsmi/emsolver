function addTestSuite(suiteName)
%addTestSuite : add a test suite in another test suite
%
%   addTestSuite(testSuiteName);
%
%
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 31/03/2005
%


if evalin('base', 'exist(''MUNIT_SUITE'', ''var'');')
    cmd = sprintf('MUNIT_SUITE.suites = {MUNIT_SUITE.suites{:}, ''%s''};', suiteName);
    evalin('base', cmd);
end
