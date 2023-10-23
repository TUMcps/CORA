function res = getVariables(obj)
% getVariables - alphabetically sorted list of variables from stl object
%
% Syntax:
%    res = getVariables(obj)
%
% Inputs:
%    obj - logic formula (class stl)
%
% Outputs:
%    res - sort list of variables
%
% Example: 
%    x = stl('x',15);
%    eq = until(x(2) < -1,x(15) < -1,interval(0,1));
%    vars = getVariables(eq);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Authors:       Niklas Kochdumper, Benedikt Seidl
% Written:       09-November-2022 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    res = []; vars = obj.variables;
    len = cellfun(@(x)length(x),vars);

    for i = 1:max(len)
       res = [res; sort(vars(len == i))];
    end
end

% ------------------------------ END OF CODE ------------------------------
