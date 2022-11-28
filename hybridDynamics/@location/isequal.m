function res = isequal(loc1,loc2)
% isequal - checks if two locations are equal by comparing the invariants,
%    transitions, flow equations, and names
%
% Syntax:  
%    res = isequal(loc1,loc2)
%
% Inputs:
%    trans1 - location object
%    trans2 - location object
%
% Outputs:
%    res - true/false
%
% Example:
%    % invariant
%    polyOpt = struct('A',[-1,0],'b',0);
%    inv = mptPolytope(polyOpt);
%    
%    % transition
%    c = [-1;0]; d = 0; C = [0,1]; D = 0;
%    guard = conHyperplane(c,d,C,D);
%
%    % reset function
%    reset = struct('A',[1,0;0,-0.75],'c',[0;0]);
%
%    % transition
%    trans{1} = transition(guard,reset,2);
%
%    % flow equation
%    dynamics = linearSys([0,1;0,0],[0;0],[0;-9.81]);
%
%    % define locations
%    loc1 = location('S1',inv,trans,dynamics);
%    loc2 = location('S2',inv,trans,dynamics);
%
%    % comparison
%    res = isequal(loc1,loc1);
%    res = isequal(loc1,loc2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      26-November-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% easy checks first to facilitate quick exits
res = true;

% compare names
if ~strcmp(loc1.name,loc2.name)
    res = false; return
end

% compare invariant
if ~isequal(loc1.invariant,loc2.invariant)
    res = false; return
end

% compare flow equations
if ~isequal(loc1.contDynamics,loc2.contDynamics)
    res = false; return
end

% compare transitions
% same number of outgoing transitions
if length(loc1.transition) ~= length(loc2.transition)
    res = false; return
end
% try to find match between transitions
idxInLoc2 = false(length(loc1.transition));
for i=1:length(loc1.transition)
    found = false;
    for j=1:length(loc2.transition)
        if ~idxInLoc2(j)
            if isequal(loc1.transition{i},loc2.transition{j})
                found = true;
                idxInLoc2(j) = true;
                break
            end
        end
    end
    if ~found
        % i-th transition in loc1 has no match in loc2
        res = false; return
    end
end

%------------- END OF CODE --------------
