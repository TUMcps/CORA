function [res,msg] = c_HA_timeStep(val,sys,options)
% c_HA_timeStep - costum validation function for options.timeSteps
%
% Syntax:
%    [res,msg] = c_HA_timeStep(val,sys,options)
%
% Inputs:
%    val - value for given param / option
%    sys - hybridAutomaton object
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
% Written:       04-February-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume check ok
res = true;
msg = '';

if ~iscell(val)
    % single value given (for all locations)
    if ~isscalar(val)
        msg = 'has to be a scalar value';
        res = false; return;
    elseif ~isnumeric(val)
        msg = 'has to be a numeric value';
        res = false; return;
    elseif val <= 0
        msg = 'must not be smaller or equal to 0';
        res = false; return;
    end
    
else
    numLoc = length(sys.location);
    % cell array, different timeStep for each location
    if ~all(size(val) == [numLoc,1]) && ~all(size(val) == [1,numLoc])
        msg = 'has to match the number of locations';
        res = false; return;
    else
        for i=1:numLoc
            if ~isscalar(val{i})
                msg = 'has to be a scalar value';
                res = false; return;
            elseif ~isnumeric(val{i})
                msg = 'has to be a numeric value';
                res = false; return;
            elseif val{i} <= 0
                msg = 'must not be smaller or equal to 0';
                res = false; return;
            end
        end
    end 
        
end


end

% ------------------------------ END OF CODE ------------------------------
