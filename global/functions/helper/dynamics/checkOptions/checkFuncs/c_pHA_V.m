function [res,msg] = c_pHA_V(val,sys,params)
% c_pHA_V - costum validation function for params.V
%
% Syntax:
%    [res,msg] = c_pHA_V(val,sys,params)
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
% Written:       31-August-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume check ok
res = true;
msg = '';

if ~iscell(val)
    % same input everywhere
    if dim(val) ~= sys.nrOfNoises
        res = false;
        msg = 'has to match the number of inputs';
        return;
    end

else
    
    % number of components in parallel hybrid automaton
    numComps = length(sys.components);
    
    if ~(all(size(val) == [numComps,1]) || all(size(val) == [1,numComps]))
        res = false;
        msg = 'has to match the number of components';
        return;
    end
    
    % go over all components
    for i=1:numComps
        numLoc = length(sys.components(i).location);
        for j=1:numLoc
            if dim(val{i}{j}) ~= sys.components(i).location(j).contDynamics.nrOfNoises
                res = false;
                msg = 'has to match number of noises for each location';
                return;
            end
        end
    end
    
end

end

% ------------------------------ END OF CODE ------------------------------
