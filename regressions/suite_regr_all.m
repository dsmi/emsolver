function suite_regr_all
% suite_regr_all : testsuite including all the solver regression tests.
%

addTest('regr_pcloop');
addTest('regr_loop');
addTest('regr_pcdipolea');
addTest('regr_dipolea');
addTest('regr_wire');
addTest('regr_tline');
addTest('regr_capacitor');
