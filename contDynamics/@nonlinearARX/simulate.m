function [tVec,ind1,ind2,y] = simulate(nlnARX,params,varargin)
% simulate - simulates a nonlinear ARX system
%
% Syntax:
%    [t,ind1] = simulate(nlnARX,params)
%    [t,ind1,ind2,y] = simulate(nlnARX,params)
%
% Inputs:
%    nlnARX - nonlinearARX object
%    params - struct containing the parameters for the simulation
%       .tStart: initial time
%       .tFinal: final time
%       .y_init: first p measurements after time tStart
%       .u: piecewise constant input signal u(t) specified as a matrix
%           for which the number of rows is identical to the number of
%           time steps (if one row: input is constant over all time
%           steps)
%
% Outputs:
%    t - time vector
%    ind1 - [] (argument exists only for syntax consistency)
%    ind2 - [] (argument exists only for syntax consistency)
%    y - output vector with number of rows = number of time steps
%
% Example: 
%    f = @(y,u) [y(1,1) + u(1,1) - y(2,1); y(3,1) + u(2,1)*cos(y(1,1)); y(5,1) + u(4,1)*sin(y(1,1))];
%    dt = 0.25;
%    nlnARX = nonlinearARX(f,dt,3,2,2)
%    
%    params.x0 = [0 1;0 1;0 1];
%    params.tFinal = 1;
%    params.u = [0.02 0.02 0.02 0.02 0.02; 1 2 3 4 5;];
%
%    [t,~,~,y] = simulate(nlnARX,params);
%
%    plot(y(:,2),y(:,3),'.k','MarkerSize',20);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSysDT/simulate

% Authors:       Laura Luetzow
% Written:       05-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set default output arguments
ind1 = [];
ind2 = [];

% parse input arguments
if ~isfield(params,'tStart')
	params.tStart = 0; 
end

tVec = (params.tStart:nlnARX.dt:params.tFinal)';

% consider changing inputs
change = false;

if size(params.u,2) ~= 1
    change = true;
    if size(params.u,2) < length(tVec)
        throw(CORAerror('CORA:specialError',...
            'Input signal "params.u" has the wrong dimension!'));
    end
end

% initialize
p = nlnARX.n_p;
if ~isfield(params, 'y0')
    params.y0 = reshape(params.x0,nlnARX.nrOfOutputs,[]);
end
y = [params.y0 zeros(nlnARX.nrOfOutputs, length(tVec)-size(params.y0,2))];

% loop over all time steps
for k = size(params.y0,2):length(tVec)-1
    x_prev = reshape(y(:,k-p+1:k),[],1);
    if change
        u_prev = reshape(params.u(:,k-p+1:k+1),[],1);
    else
        u_prev = repmat(params.u, p+1, 1);
    end

    y(:, k+1) = nlnARX.mFile(x_prev,u_prev);
end
y=y';

% ------------------------------ END OF CODE ------------------------------
