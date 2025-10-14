function traj = addNoise(traj,std,vars,varargin)
% addNoise - Add random measurement noise to the states of trajectories
%
% Syntax:
%    traj = addNoise(traj,std,vars)
%    traj = addNoise(traj,std,vars,type)
%
% Inputs:
%    traj - trajectory object
%    std - standard deviation of the noise that is added
%    vars - variables ('x', 'y', 'u' and combinations) where noise should
%               be added
%    type - type of noise that is addes ('uniform' (default) or 'gaussian')
%
% Outputs:
%    traj - trajectory object with noise
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Niklas Kochdumper, Laura Luetzow
% Written:       06-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
narginchk(3,4);
type = setDefaultValues({'uniform'},varargin);

inputArgsCheck({{traj,'att','trajectory'}, ...
    {std,'att','numeric','positive'}, ...
    {type,'str',{'uniform','gaussian'}}});

% iterate through all trajectory objects
for m = 1:length(traj)
    for var = vars
        if ~any(strcmp(var, {'x','y','u','a'}))
            throw(CORAerror('CORA:specialError','Noise can only be added to x, y, u, or a'))
        end
        % generate noise
        if strcmp(type,'uniform')
            noise = std*sqrt(3)*2*(rand(size(traj(m).(var)))-0.5);
        elseif strcmp(type,'gaussian')
            noise = normrnd(0,std,size(traj(m).(var)));
        end

        % add noise to the data
        traj(m).(var) = traj(m).(var) + noise;
        if strcmp(var, 'x')
            % initial state must be the same for all realizations
            traj(m).x(:,1,:) = repmat(traj(m).x(:,1,1), 1,1,size(traj(m).x,3));
        end
    end
end
end

% ------------------------------ END OF CODE ------------------------------
