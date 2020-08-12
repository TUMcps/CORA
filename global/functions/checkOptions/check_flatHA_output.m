function check_flatHA_output(obj)
% check_flatHA_output - check if outputs, which are not yet supported for
%                       hybrid automata, are provided
%
% Syntax:
%    options = check_flatHA_output(obj)
%
% Inputs:
%    obj     - hybrid automaton object
%
% Outputs:
%    none
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: hybridAutomaton

% Author:       Niklas Kochdumper
% Written:      15-June-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    locations = obj.location;
    
    % loop over all locations
    for i = 1:length(locations)
       sys = locations{i}.contDynamics;
       if isa(sys,'linearSys')
          if (~isempty(sys.C) && (~isscalar(sys.C) || sys.C ~= 1)) || ...
              ~isempty(sys.D) || ~isempty(sys.k)
             error('Outputs are not supported for hybrid automata!');  
          end
       end
    end
end

%------------- END OF CODE --------------