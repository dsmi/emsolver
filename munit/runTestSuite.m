function varargout = runTestSuite(testName, varargin)
%runTestSuite : run a matlab test suite
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

% check verbosity option
verbose = 1;
if length(varargin)>0
    var = varargin{1};
    if strcmp(var, 'silent') | strcmp(var, 'quiet')
        verbose = 0;
    end
end

% set MUNIT_SUITE flags and options
if ~evalin('base', 'exist(''MUNIT_SUITE'', ''var'');');
                    
    evalin('base', sprintf('MUNIT_SUITE.verbose=%d;', verbose));
    evalin('base', 'MUNIT_SUITE.passed=0;');
    evalin('base', 'MUNIT_SUITE.tested=0;');
    evalin('base', 'MUNIT_SUITE.tests={};');
    evalin('base', 'MUNIT_SUITE.suites={};');
else
    firstTest = false;
end


% run the testSuite file. Names of test to run will be stored in the field
% tests of MUNIT_SUITE
evalin('base', testName);


% if some testsuites have been added to the current test suite, process
% each test suite to find all test cases. Process until no more test suites
% are stored in variable MUNIT_SUITE.suites
suites = evalin('base', 'MUNIT_SUITE.suites;');
while length(suites)>0
    suiteName = suites{1};
    evalin('base', 'MUNIT_SUITE.suites={MUNIT_SUITE.suites{2:end}};');
    evalin('base', suiteName);
    suites = evalin('base', 'MUNIT_SUITE.suites;');
end


% get name of all the tests, and remove doubles
tests = evalin('base', 'MUNIT_SUITE.tests');
tests = unique(tests);

for i=1:length(tests)
    
    % run the test
    [p t] = eval('base', sprintf('runTest(''%s'', ''quiet'');', tests{i}));

    % update number of test passed and tested of the suite
    evalin('base', sprintf('MUNIT_SUITE.tested=MUNIT_SUITE.tested+%d;', t));
    evalin('base', sprintf('MUNIT_SUITE.passed=MUNIT_SUITE.passed+%d;', p));

    % display informations on the current test if specified
    if verbose
        if p>1
            plur = 's';
        else
            plur = ' ';
        end
        disp(sprintf('  > %-40s : %3d test%c passed on %3d', tests{i}, p, plur, t));
    end
end


% update return variables
tested = evalin('base', 'MUNIT_SUITE.tested;');
passed = evalin('base', 'MUNIT_SUITE.passed;');

% sum up result
disp(sprintf('Test suite result : %d tests passed on %d, in %d files', passed, tested, length(tests)));

% clean up global workspace
if firstTest
    evalin('base', 'clear MUNIT_SUITE;');
end    

if nargout==2
    varargout{1} = passed;
    varargout{2} = tested;
end