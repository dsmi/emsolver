function assertTrue(var)
%assertTrue : assert a variable is true (different to zero)
%
%   assertTrue(x) : give error message if x is false
%
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 01/04/2004.
%

if ~var
    % get debug information
    [stack, si] = dbstack;
    M = length(stack);

    disp(sprintf('Assertion not verified in function %s (file "%s", l-%d)', ...
    stack(2).name, stack(2).file, stack(2).line));

    % check existence of MUNIT workspace
    if evalin('base', 'exist(''MUNIT'', ''var'');')
                
        %if evalin('base', 'MUNIT.verbose');
        %end
        
        % set the assert flag to false
        evalin('base', 'MUNIT.assertFlag=0;');
    end
end