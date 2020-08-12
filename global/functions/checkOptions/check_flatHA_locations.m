function whichCheck = check_flatHA_locations(options, obj)
% check_flatHA_locations - checks if all the locations have the necessary
%   options for their checkOptionsReach
%
% Syntax:
%    check_flatHA_locations(options, obj)
%
% Inputs:
%    options   - options for object
%    obj       - system (e.g. linearSys)
%
% Outputs:
%    whichCheck - cells of strings which systems are in locations
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
% Written:      06-Mar-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

locations = obj.location;
numLoc = length(locations);


% depending on the locations inside the hybridAutomaton object,
%  their options have to be tested:
%    - linearSys
%    - linParamSys
%    - nonlinDASys
%    - nonlinearSys
%    - nonlinParamSys
systems = {'linearSys'; 'linParamSys';...
    'nonlinDASys'; 'nonlinearSys'; 'nonlinParamSys'};
% define whichCheck: true / false for everyone of the below listed systems
whichCheck = false(1,length(systems));
% save indices of locations corresponding to systems for checkOptions later
idxLocation = zeros(1,length(systems));

% determine which check(s) necessary
for i=1:length(systems)
    % search along systems, if found, change corresponding index
    %  in whichCheck to 1 and search for next system
    for j=1:numLoc
        if isa(locations{j}.contDynamics, systems{i})
            whichCheck(i) = true;
            idxLocation(i) = j;
            break;
        end
    end
end

% execute checkOptionsReach for found systems in flatHA
for i=1:length(idxLocation)
    if idxLocation(i)
        % idxLocation nonzero -> corresponding system present in flatHA
        % put hyb param of checkOptionsReach to 1, as called from HA
        checkOptionsReach(locations{idxLocation(i)}.contDynamics,...
            options,1);
    end
end

end

%------------- END OF CODE --------------
