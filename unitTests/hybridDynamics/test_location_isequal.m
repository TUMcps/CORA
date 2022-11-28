function res = test_location_isequal
% test_location_isequal - test function for isequal
%
% Syntax:  
%    res = test_location_isequal
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
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

% assume true
res = true;

% invariant
polyOpt = struct('A',[-1,0],'b',0);
inv1 = mptPolytope(polyOpt);
polyOpt = struct('A',[-1,0],'b',0.5);
inv2 = mptPolytope(polyOpt);

% transition
guard = conHyperplane([-1;0],0,[0,1],0);

% reset function
reset1 = struct('A',[1,0;0,-0.75],'c',[0;0]);
reset2 = struct('A',[1,0;0,-0.75],'c',[1;0]);

% transition
trans1 = transition(guard,reset1,1);
trans2 = transition(guard,reset1,2);
trans3 = transition(guard,reset2,2);

% flow equation
dynamics1 = linearSys([0,1;0,0],[0;0],[0;-9.81]);
dynamics2 = linearSys([0,1;0,0],[0;0],[0;-9.80]);


% same location
if ~isequal(location(inv1,{trans1},dynamics1),...
            location(inv1,{trans1},dynamics1))
    res = false;
end
if ~isequal(location(inv1,{trans1,trans2},dynamics1),...
            location(inv1,{trans2,trans1},dynamics1))
    res = false;
end

% different invariant
if isequal(location(inv1,{trans1},dynamics1),...
           location(inv2,{trans1},dynamics1))
    res = false;
end

% different dynamics
if isequal(location(inv1,{trans1},dynamics1),...
           location(inv1,{trans1},dynamics2))
    res = false;
end

% different transition
if isequal(location(inv1,{trans1},dynamics1),...
           location(inv1,{trans2},dynamics1))
    res = false;
end

% different transition
if isequal(location(inv1,{trans1},dynamics1),...
           location(inv1,{trans1,trans3},dynamics1))
    res = false;
end

%------------- END OF CODE --------------
