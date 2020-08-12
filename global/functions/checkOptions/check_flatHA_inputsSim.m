function options = check_flatHA_inputsSim(options, obj)
% check_flatHA_inputsSim - checks if inputs (u)
%  1) exist
%  2) take an allowed value
%
% Syntax:
%    options = check_flatHA_inputsSim(options, obj)
%
% Inputs:
%    options - options for object
%    obj     - hybrid automaton object
%
% Outputs:
%    options   - options for object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   -

% Author:       Mark Wetzlinger
% Written:      13-Mar-2019
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

locations = obj.location;
numLoc = length(locations);

% check if Uloc is specified -> warning that it is redundant
if isfield(options,'uLoc')
    warning('options.uLoc is redundant!');
end

option = 'u';
% assign internal value options.Uloc which stores the input set for each
% location
if ~isfield(options,option)
% initialize input set with zeros

    uloc = cell(numLoc,1);
    
    for i = 1:numLoc
        sys = locations{i}.contDynamics;
        uloc{i} = zeros(max(1,sys.nrOfInputs),1);
    end    
    
elseif ~iscell(options.u) 
% same input set for each location    
    
    uloc = cell(numLoc,1);

    for i = 1:numLoc
        uloc{i} = options.u;
    end
    
    options = rmfield(options,option);
    
else
% copy input set for each location

    if all(size(options.u) ~= [numLoc,1]) && ...
       all(size(options.u) ~= [1,numLoc])
        error(printOptionOutOfRange(obj,option,strct));
    else
        uloc = options.u;
    end
    
    options = rmfield(options,option);
end

options.uLoc = uloc;

end
%------------- END OF CODE --------------