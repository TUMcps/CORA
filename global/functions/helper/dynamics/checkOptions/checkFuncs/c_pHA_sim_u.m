function [res,msg] = c_pHA_sim_u(val,sys,params)
% c_pHA_sim_u - costum validation function for params.u
%
% Syntax:
%    [res,msg] = c_pHA_sim_u(val,sys,options)
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

if ~iscell(val)
    % same input everywhere
    if length(val) ~= sys.nrOfInputs
        msg = 'has to match the number of inputs';
        res = false; return;
    end

else
    
    numComps = length(sys.components);
    
    if ~(all(size(val) == [numComps,1]) || all(size(val) == [1,numComps]))
        msg = 'has to match the number of components';
        res = false; return;
    end
    
    % go over all components
    for i=1:numComps
    	% check if size of input matches inputCompMap
    	dims = length(find(params.inputCompMap == i));
        if dims > 0
            numLoc = length(sys.components(i).location);
            if ~(all(size(val{i}) == [numLoc,1]) || all(size(val{i}) ~= [1,numLoc]))
                msg = 'has to match the number of locations';
                res = false; return;
            end
            for j=1:numLoc
                if size(val{i}{j},1) ~= dims
                    msg = 'has to match inputCompMap';
                    res = false; return;
                end
            end
        end
    end
    
end

end

% ------------------------------ END OF CODE ------------------------------
