function [t,ind1,ind2,y] = simulate(linARX,params,varargin)
% simulate - simulates a linear ARX system
%
% Syntax:
%    [t,ind1,ind2,y] = simulate(linARX,params)
%
% Inputs:
%    linARX - linearARX system
%    params - struct containing the parameters for the simulation
%       .tStart: initial time 
%       .tFinal: final time
%       .y0: first p measurements after time tStart
%       .u: piecewise-constant input signal u(t) of dimension
%               dim_u x z, where
%               z = 1 (fixed input) or
%                   N (if B_bar{1}=0)
%                   N+1 
%                   where N is the number of steps
%
% Outputs:
%    t - time vector
%    ind1 - [] (argument exists only for syntax consistency)
%    ind2 - [] (argument exists only for syntax consistency)
%    y - output vector
%
% Example:
%    dt = 0.1;
%    A_bar = {[-0.4 0.6; 0.6 -0.4];[0.1 0; 0.2 -0.5]};
%    B_bar = {[0; 0];[0.3; -0.7];[0.1; 0]};
%    linARX = linearARX(A_bar,B_bar,dt)
%
%    params.y0 = [[0;0] [0.1;0.2]];
%    params.tFinal = 0.4;
%    params.u = [0.1 0.05 0.05 -0.1];
%
%    [t,~,~,y] = simulate(linARX,params);
%
%    plot(y(1,:),y(2,:),'.k','MarkerSize',20);
%
% Reference:
%   [1] L. Luetzow, M. Althoff. "Reachability Analysis of ARMAX Models", in
%       Proc. of the 62th International Conference on Decision and Control,
%       2023.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Laura Luetzow
% Written:       02-February-2023
% Last update:   28-August-2025 (LL, transpose t and y)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
if ~isfield(params,'tStart')
	params.tStart = 0; 
end

% set default output arguments
ind1 = [];
ind2 = [];

% compute time vector and number of steps
p = linARX.n_p;
t = params.tStart:linARX.dt:params.tFinal;

% check u, w, and v; extend fixed vector to matrices for loop below
params = aux_uwv(linARX,params,length(t)-1);

% initialize
if ~isfield(params, 'y0')
    params.y0 = reshape(params.x0,linARX.nrOfOutputs,[]);
end
y = [params.y0 zeros(size(params.y0,1), length(t)-p)];

if ~linARX.tvp % simulate outputs with normal ARX parametrization
    % loop over all time steps
    for k=p:length(t)-1
        y(:, k+1) = linARX.B_bar{1}*params.u(:, k+1);
        for i=1:p
            y(:, k+1) = y(:, k+1) + linARX.A_bar{i,1} * y(:, k+1-i) ...
                + linARX.B_bar{i+1} * params.u(:, k+1-i);
        end
    end
else % simulate outputs with time-varying parametrization
    y_init = reshape(params.y0, [],1);
    k=p;
    while k<=length(t)-1
        k_plus = k + p - 1;
        if ~isfield(linARX, "A_tilde") || length(linARX.A_tilde) < k+1 || isempty(linARX.A_tilde(k+1))
            % copmute tyime-varying parameters
            computeTVP(linARX,k, 0:k_plus);
        end

        % compute output at time step k
        y_add = linARX.A_tilde{k+1} * y_init;
        for i = 0:k_plus
            if k_plus -i <= length(t)-1
                y_add = y_add + linARX.B_tilde{i+1,k+1} * params.u(:,k_plus + 1 -i);
            end
        end
        y(:, k+1:k_plus+1) = reshape(y_add,linARX.nrOfOutputs,[]);

        % increment time step k by p
        k = k+p;
    end
    y = y(:, 1:length(t));
end
end


% Auxiliary functions -----------------------------------------------------

function params = aux_uwv(linARX,params,steps)

% check input, set default
if ~isfield(params,'u')
    params.u = zeros(linARX.nrOfInputs,steps+1);
else
    % length always steps+1
    if size(params.u,2) == 1
        params.u = repmat(params.u,1,steps+1);
    elseif size(params.u,2) == steps && ~any(any(linARX.B_bar{1}))
        params.u = [params.u, zeros(size(params.u,1),1)];
    elseif size(params.u,2) ~= steps+1
        throw(CORAerror('CORA:wrongFieldValue','params.u',...
            'has to be of size(2) of steps+1'));
    end
end
end

% ------------------------------ END OF CODE ------------------------------
