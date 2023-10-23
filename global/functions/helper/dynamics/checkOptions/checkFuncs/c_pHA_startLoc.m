function [res,msg] = c_pHA_startLoc(val,sys,params)
% c_pHA_startLoc - costum validation function for params.startLoc
%
% Syntax:
%    [res,msg] = c_pHA_startLoc(val,sys,options)
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

% Authors:       Mark Wetzlinger
% Written:       04-February-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume check ok
res = true;
msg = '';

numComp = length(sys.components);

if ~all(size(params.startLoc) == [numComp,1])
    msg = 'has to match the number of components';
    res = false; return;
else
    for c=1:numComp
        % loop through every component
        comp = sys.components(c);
        numLoc = length(comp.location);
        if params.startLoc(c) > numLoc
            msg = 'must not be larger than the number of locations';
            res = false; return;
        elseif params.startLoc(c) <= 0
            msg = 'must not be smaller or equal to 0';
            res = false; return;
        elseif mod(params.startLoc(c),1.0) ~= 0
            msg = 'has to be an integer value';
            res = false; return;
        end
    end
end


end

% ------------------------------ END OF CODE ------------------------------
