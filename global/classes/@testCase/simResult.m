function simRes = simResult(obj)
% simResult - converts a test case to a simResult object.
%
% Syntax:
%    simRes = simResult(obj)
%
% Inputs:
%    obj - testCase object
%
% Outputs:
%    simRes - simResult object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: @simResult

% Authors:       Matthias Althoff
% Written:       21-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% create time vector
% number of time steps
nrOfTimeSteps = size(obj.y,1);
% time incerement
dt = obj.sampleTime;
% time vector
tVec = 0:dt:dt*(nrOfTimeSteps-1);

% create cell arrays with one entry only
x{1} = obj.x;
t{1} = tVec';
y{1} = obj.y;

% create simResult object
simRes = simResult(x,t,{},y);

% ------------------------------ END OF CODE ------------------------------
