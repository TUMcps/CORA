function res = isequal(loc1,loc2,varargin)
% isequal - checks if two locations are equal by comparing the invariants,
%    transitions, flow equations, and names
%
% Syntax:
%    res = isequal(loc1,loc2)
%    res = isequal(loc1,loc2,tol)
%
% Inputs:
%    trans1 - location object
%    trans2 - location object
%    tol - tolerance (optional)
%
% Outputs:
%    res - true/false
%
% Example:
%    % invariant
%    inv = polytope([-1,0],0);
%    
%    % transition
%    c = [-1;0]; d = 0; C = [0,1]; D = 0;
%    guard = conHyperplane(c,d,C,D);
%
%    % reset function
%    reset = struct('A',[1,0;0,-0.75],'c',[0;0]);
%
%    % transition
%    trans = transition(guard,reset,2);
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

% Authors:       Mark Wetzlinger
% Written:       26-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% too many input arguments
if nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% default values
tol = setDefaultValues({eps},varargin);

% check input arguments
inputArgsCheck({{loc1,'att','location'};
                {loc2,'att','location'};
                {tol,'att','numeric',{'scalar','nonnegative','nonnan'}}});

% check array length
if any(size(loc1) ~= size(loc2))
    res = false; return
end

% loop over all individual transition objects
[r,c] = size(loc1);
res = false(r,c);
for i=1:r
    for j=1:c
        res(i,j) = aux_isequal(loc1(i,j),loc2(i,j),tol);
    end
end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isequal(loc1,loc2,tol)

% assume true
res = true;

% compare names
% if ~strcmp(loc1.name,loc2.name)
%     res = false; return
% end

% compare invariants: 
% note: we use 'eq' instead of 'isequal' as long as the polytope class
% exists, once the switch to the polytope class is done, use 'isequal' in
% accordance with other calls
if any([isnumeric(loc1.invariant),isnumeric(loc2.invariant)])
    % empty location object may have .invariant = []
    if xor(isnumeric(loc1.invariant),isnumeric(loc2.invariant))
        res = false; return
    end
elseif ~eq(loc1.invariant,loc2.invariant,tol)
    res = false; return
end

% compare flow equations
if ~isequal(loc1.contDynamics,loc2.contDynamics,tol)
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
    % assume no matching transition was found
    found = false;

    % loop over all transitions of second location
    for j=1:length(loc2.transition)
        % skip transitions that have already been matched
        if ~idxInLoc2(j)
            % check for equality
            if isequal(loc1.transition(i),loc2.transition(j),tol)
                % matching transition found
                found = true; idxInLoc2(j) = true;
                break
            end
        end
    end
    
    if ~found
        % i-th transition in loc1 has no match in loc2
        res = false; return
    end
end

end

% ------------------------------ END OF CODE ------------------------------
