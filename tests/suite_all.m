function suite_all
% suite_all : testsuite including all the solver components tests.
%

addTestSuite('suite_intg');
addTestSuite('suite_common');
addTestSuite('suite_mesh');
addTestSuite('suite_solver');
addTestSuite('suite_tlines');
addTestSuite('suite_layers');
addTestSuite('suite_mom2d');
