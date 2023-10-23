function [res,msg] = c_pHA_R0(val,sys,params)
% c_pHA_R0 - costum validation function for params.R0
%
% Syntax:
%    [res,msg] = c_pHA_R0(val,sys,params)
%
% Inputs:
%    val - value for given param / option
%    sys - parallelHybridAutomaton object
%    params - model parameters
%
% Outputs:
%    res - logical whether validation was successful
%    msg - error message if validation failed
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   -

% Authors:       Mark Wetzlinger, Niklas Kochdumper
% Written:       04-February-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume check ok
res = true;
msg = '';

totalStates = 0;
for i=1:length(sys.bindsStates)
    totalStates = totalStates + length(sys.bindsStates{i});
end
if dim(params.R0) ~= totalStates
    msg = 'has to match the number of state binds';
    res = false;
end

end

% ------------------------------ END OF CODE ------------------------------
