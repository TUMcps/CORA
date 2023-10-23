function [dt,uniform,hybrid] = timeStepSize(R)
% timeStepSize - returns time step size for the reachable set
%
% Syntax:
%    [dt,uniform,hybrid] = timeStepSize(R)
%
% Inputs:
%    R - reachSet object
%
% Outputs:
%    dt - time step size (scalar or vector)
%    uniform - uniform time step size (true) or not (false)
%    hybrid - reachable set belongs to a hybrid system
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachSet

% Authors:       Niklas Kochdumper
% Written:       08-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input argument
inputArgsCheck({{R,'att',{'reachSet'},{''}}});

% initialization
hybrid = false;
uniform = false;
dt = [];

% loop over all reachable set objects
for i = 1:size(R,1)

    % check if times are intervals or scalar values
    num = cellfun(@isnumeric,R(i,1).timePoint.time);

    if all(num)
        dt = [dt; diff(cell2mat(R(i,1).timePoint.time))];
    else
        tmp = cellfun(@infimum,R(i,1).timePoint.time);
        dt = [dt; diff(tmp)];
        hybrid = true;
    end
end

% check if time step is uniform
if all(abs(diff(dt)) < 1e-10)
    uniform = true;
    dt = dt(1);
end

% ------------------------------ END OF CODE ------------------------------
