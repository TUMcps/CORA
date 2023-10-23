function [res,msg] = c_tu(val,sys,params,options)
% c_tu - costum validation function to check params.tu
%
% Syntax:
%    [res,msg] = c_tu(val,sys,params,options)
%
% Inputs:
%    val - value for options.timeStep
%    sys - some contDynamics object
%    params - model parameters
%    options - algorithm parameters
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

% Authors:       Mark Wetzlinger
% Written:       25-May-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume check ok
res = true;
msg = '';

if val == 0
    % default value
    return;
% elseif size(val,2) ~= size(params.u,2)
%     % check if size of params.tu matches params.u
%     res = false;
%     msg = ['has to match the length of params.u'];
%     return;
elseif val(1) ~= params.tStart
    % first entry has to be start time
    res = false;
    msg = ['''s first entry has to match params.tStart'];
    return;
elseif ~withinTol(val(end),params.tFinal) && val(end) > params.tFinal
    % last entry cannot be larger than final time
    res = false;
    msg = ['''s last entry has to be smaller than params.tFinal'];
    return;
end

end

% ------------------------------ END OF CODE ------------------------------
