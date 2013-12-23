function varargout = runTest(testName, varargin)
%runTest : run a matlab unit test
%
%   usage :
%   [p, t] = runTest('myTest');
%   lauch the program 'myTest', count the number of assertions non
%   verified, and return [0 1] or [1 1] whether all assertion are verified
%   or not.
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


firstTest = true;
verbose = true;
if length(varargin)>0
    var = varargin{1};
    if strcmp(var, 'silent') | strcmp(var, 'quiet')
        verbose = false;
    end
end

% set MUNIT flags if first test
if ~evalin('base', 'exist(''MUNIT'', ''var'');');
    evalin('base', sprintf('MUNIT.verbose=%d;', verbose));
    evalin('base', 'MUNIT.passed=0;');
    evalin('base', 'MUNIT.tested=0;');
else
    firstTest=false;
end


% set flag to true
evalin('base', 'MUNIT.assertFlag=1;');

% launch the test. The test will call 'assertXXX' functions, which will
% possibly set MUNIT_ASSERT_FLAG to false.
ind = strfind(testName, '/');
if length(ind)>0
    base = pwd;
    path = testName(1:ind(end)-1);
    name = testName(ind(end)+1:end);
    
    cd(path);
    eval(name);
    cd(base);
else
    eval(testName);
end


% update number of passed and run tests
evalin('base', 'MUNIT.tested=MUNIT.tested+1;');
if evalin('base', 'MUNIT.assertFlag');
    evalin('base', 'MUNIT.passed=MUNIT.passed+1;');
end


% update return variables
tested = evalin('base', 'MUNIT.tested;');
passed = evalin('base', 'MUNIT.passed;');

% show result if specified
if verbose
    disp(sprintf('Result : %d test passed on %d', passed, tested));
end


% clean up global workspace
if firstTest
    evalin('base', 'clear MUNIT;');
end    

if nargout==2
    varargout{1} = passed;
    varargout{2} = tested;
end