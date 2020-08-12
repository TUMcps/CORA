function options = check_flatHA_inputs(options,obj)
% check_flatHA_inputs - checks if inputs (U)
%  1) exist
%  2) take an allowed value
%
% Syntax:
%    options = check_flatHA_inputs(options, obj)
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

% Author:       Mark Wetzlinger, Niklas Kochdumper
% Written:      13-Mar-2019
% Last update:  03-May-2020 (change Uloc to U)
%               03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

locations = obj.location;
numLoc = length(locations);
strct = 'params';

% check if Uloc is specified -> warning that it is redundant
if isfield(options,'Uloc')
    warning('options.Uloc is redundant!');
end

option = 'U';
% assign internal value options.Uloc which stores the input set for each
% location
if ~isfield(options,option)
% initialize input set with zeros

    Uloc = cell(numLoc,1);
    
    for i = 1:numLoc
        sys = locations{i}.contDynamics;
        Uloc{i} = zonotope(zeros(max(1,sys.nrOfInputs),1));
    end    
    
elseif ~iscell(options.U) 
% same input set for each location    
    
    Uloc = cell(numLoc,1);

    for i = 1:numLoc
        Uloc{i} = options.U;
    end
    
    options = rmfield(options,option);
    
else
% copy input set for each location

    if all(size(options.U) ~= [numLoc,1]) && ...
       all(size(options.U) ~= [1,numLoc])
        error(printOptionOutOfRange(obj, option, strct));
    else
        Uloc = options.U;
    end
    
    options = rmfield(options,option);
end

options.Uloc = Uloc;

end

%------------- END OF CODE --------------

