function [guards,setIndices] = potInt(obj,R,options)
% potInt - determines which reachable sets potentially intersect with guard
%          sets
%
% Syntax:  
%    [guards,setIndices] = potInt(obj,R,options)
%
% Inputs:
%    obj - location object
%    R - cell-array of reachable sets
%
% Outputs:
%    guards - guards that are potentially intersected
%    setIndices - indices of the reachable sets that intersect the guards
%    options - struct containing the algorithm settings
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      08-May-2007 
% Last update:  26-October-2007
%               20-October-2010
%               27-July-2016
%               23-November-2017
%               03-December-2019 (NK, use approximate intersection test)
% Last revision:---

%------------- BEGIN CODE --------------

    % initialization
    N = length(R);
    M = length(obj.transition);

    guards = zeros(M*N,1);
    setIndices = zeros(M*N,1);
    
    counter = 1;
    
    % loop over all transitions
    for i = 1:M
       
        guardSet = obj.transition{i}.guard;
        target = obj.transition{i}.target;
        
        % check if terminal location is reached
        if ~all(target == options.finalLoc)
        
            % loop over all reachable sets
            for j = 1:N

                % check if reachable set intersects the guard set
                if isIntersecting(guardSet,R{j},'approx')
                    guards(counter) = i;
                    setIndices(counter) = j;
                    counter = counter + 1;
                end
            end
        end
    end
    
    % truncate the resuling lists
    guards = guards(1:counter-1);
    setIndices = setIndices(1:counter-1);

%------------- END OF CODE --------------