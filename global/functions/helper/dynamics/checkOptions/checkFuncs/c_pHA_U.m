function [res,msg] = c_pHA_U(val,sys,params)
% c_pHA_U - costum validation function for params.U
%
% Syntax:
%    [res,msg] = c_pHA_U(val,sys,options)
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
    if dim(val) ~= sys.nrOfInputs
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
    	% check if size of input matches inputCompMap
    	dims = length(find(params.inputCompMap == i));
        if dims > 0
            numLoc = length(sys.components(i).location);
            if ~(all(size(val{i}) == [numLoc,1]) || all(size(val{i}) ~= [1,numLoc]))
                res = false;
                msg = 'has to match the number of locations';
                return;
            end
            for j=1:numLoc
                if dim(val{i}{j}) ~= dims
                    res = false;
                    msg = 'has to match inputCompMap';
                    return;
                end
            end
        end
    end
    
end

end

% ------------------------------ END OF CODE ------------------------------
