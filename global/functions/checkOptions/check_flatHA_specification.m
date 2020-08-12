function options = check_flatHA_specification(options,obj,spec)
% check_flatHA_specification - brings specifications to the correct format
%
% Syntax:
%    options = check_flatHA_specification(options,obj)
%
% Inputs:
%    options - options for object
%    obj     - hybrid automaton object
%    spec    - object of class specification
%
% Outputs:
%    options   - options for object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: hybridAutomaton/reach

% Author:       Niklas Kochdumper
% Written:      07-June-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    locations = obj.location;
    numLoc = length(locations);

    % initialize specifications with empty cells
    if isempty(spec)

        specLoc = cell(numLoc,1);

    % same specification for each location  
    elseif ~iscell(spec)   

        specLoc = cell(numLoc,1);

        for i = 1:numLoc
            specLoc{i} = spec;
        end

    % copy input set for each location
    else

        if all(size(spec) ~= [numLoc,1])
            error('Input argument "spec" has the wrong format!');
        else
            specLoc = spec;
        end
    end

    options.specificationLoc = specLoc;

end

%------------- END OF CODE --------------

