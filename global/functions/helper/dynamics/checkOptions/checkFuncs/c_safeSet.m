function [res,msg] = c_safeSet(val,sys,params,options)
% c_safeSet - costum validation function to check params.safeSet;
%    this function is required because the values are structs
%
% Syntax:
%    [res,msg] = c_safeSet(val,sys,params,options)
%
% Inputs:
%    val - value for params.safeSet
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
% Written:       07-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume check ok
res = true;
msg = '';

allowedClasses = getMembers('safeSet');
timeHorizon = interval(params.tStart,params.tFinal);

% loop over all entries
for i=1:length(val)
    % check if correct set representation
    if ~any(ismember(allowedClasses,class(val{i}.set)))
        res = false;
        msg = getErrorMessage('membersafeSet');
        return
    end

    % check if dimension is correct
    if dim(val{i}.set) ~= sys.nrOfOutputs
        res = false;
        msg = getErrorMessage('eqoutput');
        return
    end

    % check if time is plausible
    if ~representsa(val{i}.time,'emptySet') && ~contains(timeHorizon,val{i}.time)
        res = false;
        msg = ' : time of set is outside of time horizon';
        return
    end
end

end

% ------------------------------ END OF CODE ------------------------------
