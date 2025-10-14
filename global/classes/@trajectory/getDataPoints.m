function points = getDataPoints(traj, compute_dx, phi)
% getDataPoints - transform the trajectory into a list of data points
%
% Syntax:
%    points = getDataPoints(traj)
%    points = getDataPoints(traj, compute_dx)
%    points = getDataPoints(traj, compute_dx, phi)
%
% Inputs:
%    traj - trajectory object
%    compute_dx - boolean determining if dx is computed
%    phi - template functions phi(x,u) for the dynamics \dot x = A*phi(x,u)
%          (function handle)
%
% Outputs:
%    points - struct containing the data points
%       .x - states
%       .xNext - next states
%       .dx - state derivatives (if compute_dx is true)
%       .u - inputs
%       .dt - time steps
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: trajectory

% Authors:       Niklas Kochdumper, Laura Luetzow
% Written:       28-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin < 2
    compute_dx = true;
end

points.x = [];
points.dx = [];
points.u = [];
points.xNext = [];
points.dt = [];

for i = 1:length(traj)
    if size(traj(i).x,3) > 1
        throw(CORAerror("CORA:notSupported", ['Trajectory with n_s>1 is ' ...
            'not implemented for getDataPoints.']))
    end

    if compute_dx
        % compute the derivatives
        dx = aux_derivative(traj(i).t,traj(i).x);
        points.dx = [points.dx dx(:,1:end-1)];
    end

    m = size(traj(i).x,2)-1;

    points.x = [points.x traj(i).x(:,1:end-1)];
    points.xNext = [points.xNext traj(i).x(:,2:end)];
    points.dt = [points.dt diff(traj(i).t)];

    if ~isempty(traj(i).u)
        points.u = [points.u traj(i).u(:,1:m)];
    end
end
if nargin == 3
    % transform data points by the template function for the dynamics
    points = aux_transformBasisFunc(points,phi);
end

end


% Auxiliary functions -----------------------------------------------------

function [dy,x] = aux_derivative(x,y)
% numerically estimate the time derivative for the data
    
    % number of subintervals
    N = length(x)-1;
    
    % preallocates vector to store derivative
    dy = zeros(size(y));
    
    % approximates derivative at lower bound using forward difference
    dy(:,1) = (y(:,2)-y(:,1))/(x(2)-x(1));
    
    % approximates derivative at upper bound using backward difference
    dy(:,N+1) = (y(:,N+1)-y(:,N))/(x(N+1)-x(N));
    
    % approximates derivatives at all other nodes using central differences
    for i = 2:N
        dy(:,i) = (y(:,i+1)-y(:,i-1))/(x(i+1)-x(i-1));
    end  
end

function points = aux_transformBasisFunc(points,phi)
% transform data points by the nonlinear template function for the dynamics

    points.phi = [];
        
    for i = 1:size(points.x,2)
        if ~isempty(points.u)
            points.phi = [points.phi,phi(points.x(:,i),points.u(:,i))];
        else
            points.phi = [points.phi,phi(points.x(:,i),0)];
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
