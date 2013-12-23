function assertEquals(var1, var2, varargin)
%assertEquals : assert equality of two variables
%
%   assertEquals(v1, v2) : assert strict equality
%
%   assertEquals(v1, v2, tol) : specify tolerance
%
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 01/04/2004.
%

tol = 0;
if length(varargin)>0
    tol = varargin{1};
end        

dif = abs(var1-var2)>tol;    
if sum(dif(:))>0
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