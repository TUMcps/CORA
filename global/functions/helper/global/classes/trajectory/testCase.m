function traj = testCase(varargin)
% testCase - (deprecated) object constructor for testCase
%
% Syntax:
%    traj = testCase(y,u,x,dt)
%    traj = testCase(y,u,x0,dt)
%    traj = testCase(y,u,x,dt,name)
%    traj = testCase(y,u,x0,dt,name)
%
% Inputs:
%    y - (a x q x s) vector of the measured outputs samples
%    u - (a x p x s) vector of input samples
%    x - (a x n x s) vector of state samples
%    x0 - (n x 1 x s) vector of initial states
%    dt - sampling time
%    name - name of the test case
%
% Outputs:
%    traj - trajectory object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: trajectory

% Authors:       Laura Luetzow
% Written:       26-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

CORAwarning("CORA:deprecated", "class", "testCase", ...
    "CORA v2026", ...
    "Please use trajectory instead. The object is returned as a trajectory.",...
    "The new trajectory class unifies the classes simResult and testCase.")

% avoid empty instantiation
if nargin == 0
    % empty object
    return
end
assertNarginConstructor(4:5,nargin);

if nargin == 1 && isa(varargin{1},'testCase')
    % convert to trajectory
    traj = testCase2traj(varargin{1});
    return
end
y = []; u = []; x = []; dt = []; name = "";
[y,u,x,dt,name] = setDefaultValues({y,u,x,dt,name},varargin);

% set initial state
if size(x,2) == 1 
    x0 = x;
    x = [];
else
    x0 = permute(x(1,:,:), [2 1 3]); % first row is initial state
end

% create dummy testCase object and convert to trajectory
testC = struct('y', {y}, 'u', {u}, 'x', {x}, 'initialState', {x0}, 'sampleTime', dt, 'name', name);

% create trajectory object
traj = testCase2traj(testC);

end

% ------------------------------ END OF CODE ------------------------------
