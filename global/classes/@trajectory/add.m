function traj = add(traj1,traj2)
% add - joins two trajectory objects, which have the same initial state,
%   same input vectors, same time vectors, and same locations vectors
%
% Syntax:
%    traj = add(traj1,traj2)
%
% Inputs:
%    traj1 - trajectory object
%    traj2 - trajectory object
%
% Outputs:
%    traj - resulting trajectory object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachSet

% Authors:       Laura Luetzow
% Written:       26-July-2025             
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% special case
if isempty(traj1)
     traj = traj2;
     return;
elseif isempty(traj2)
     traj = traj1;
     return;
end

% check input arguments
inputArgsCheck({{traj1,'att',{'trajectory'}};
                {traj2,'att',{'trajectory'}}});

% general case
if any(traj1.loc ~= traj2.loc, 'all') || ... % different locations
        traj1.n_k ~= traj2.n_k || ... % different time horizon
        any(traj1.t ~= traj2.t, 'all') || ... % different time steps
        any(traj1.u ~= traj2.u, 'all') || ... % different inputs
        (~isempty(traj1.x) && ~isempty(traj2.x) && any(traj1.x(:,1,1) ~= traj2.x(:,1,1), 'all')) || ... % different initial states
        any(size(traj1.x,1,2,3) ~= size(traj2.x,1,2,3), 'all') || ... % different state dimensions
        any(size(traj1.y,1,2,3) ~= size(traj2.y,1,2,3), 'all') || ... % different output dimension
        any(size(traj1.a,1,2,3) ~= size(traj2.a,1,2,3), 'all') % different algebraic dimension
    throw(CORAerror('CORA:specialError','Objects are not compatible!')); 
end

% stack states
if size(traj1.x, 2) == 1
    % only initial state is known
    x = traj1.x;
else
    % observed states are concatenated
    x = cat(3, traj1.x, traj2.x);
end

% same inputs
u = traj1.u;

% stack outputs
y = cat(3, traj1.y, traj2.y);

% same time vector
t = traj1.t;

% stack algebraic variables
a = cat(3, traj1.a, traj2.a);

% same locations
loc = traj1.loc;

% combine names
name = sprintf('%s %s', traj1.name, traj2.name);

traj = trajectory(u,x,y,t,[],loc,a,name);

% ------------------------------ END OF CODE ------------------------------
