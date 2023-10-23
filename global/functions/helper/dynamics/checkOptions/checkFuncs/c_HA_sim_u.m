function [res,msg] = c_HA_sim_u(val,sys,params)
% c_HA_sim_u - custom validation function for params.U
%
% Syntax:
%    [res,msg] = c_HA_sim_u(val,sys,options)
%
% Inputs:
%    val - value for given param / option
%    sys - hybridAutomaton object
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

if ~iscell(params.u)
    
    % check if input has the correct format
    if size(params.u,2) ~= 1
        msg = 'has to be a vector';
        res = false;
    end
	
else
    % input is a cell
    
    locations = sys.location;
    numLoc = length(locations);
    % check if input matches the number of locations
    if all(size(params.u) ~= [numLoc,1]) && ...
       all(size(params.u) ~= [1,numLoc])
        msg = 'has to match the number of locations';
        res = false; return;
    end
    
    % check input for each location
    for i = 1:length(params.u)
        if size(params.u{i},2) ~= 1
            msg = 'has to be a vector';
            res = false; return;
        end
    end
    
end


end

% ------------------------------ END OF CODE ------------------------------
