function [res,msg] = c_pHA_finalLoc(val,sys,params)
% c_pHA_finalLoc - costum validation function for params.finalLoc
%
% Syntax:
%    [res,msg] = c_pHA_finalLoc(val,sys,options)
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

% number of components
numComp = length(sys.components);

if ~all(size(params.finalLoc) == [numComp,1])
    msg = 'has to match the number of locations';
    res = false; return;
else
    for c=1:numComp
        % loop through every component
        comp = sys.components(c);
        numLoc = length(comp.location);
        if params.finalLoc(c) > numLoc + 1
            msg = 'must not be larger than the number of locations plus one';
            res = false; return;
        elseif params.finalLoc(c) < 0
            msg = 'must not be smaller or equal to 0';
            res = false; return;
        elseif mod(params.finalLoc(c),1.0) ~= 0
            msg = 'has to be an integer value';
            res = false; return;
        end
    end
end


end

% ------------------------------ END OF CODE ------------------------------
