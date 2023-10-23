function [t,ind1,ind2,y] = simulate(sys,params)
% simulate - simulates a linear ARMAX system
%
% Syntax:
%    [t,ind1,ind2,y] = simulate(sys,params)
%
% Inputs:
%    sys - linearARMAX system
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
%    sys = linearARMAX(A_bar,B_bar,dt)
%
%    params.y0 = [[0;0] [0.1;0.2]];
%    params.tFinal = 0.4;
%    params.u = [0.1 0.05 0.05 -0.1];
%
%    [t,~,~,y] = simulate(sys,params);
%
%    plot(y(:,1),y(:,2),'.k','MarkerSize',20);
%
% Reference:
%   [1] L. Luetzow, M. Althoff. "Reachability Analysis of ARMAX Models", in
%       Proc. of the 62th International Conference on Decision and Control,
%       2023.

% Authors:       Laura Luetzow
% Written:       02-February-2023
% Last update:   ---
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
p = sys.dim;
t = (params.tStart:sys.dt:params.tFinal)';

% check u, w, and v; extend fixed vector to matrices for loop below
params = aux_uwv(sys,params,length(t)-1);

% initialize
y = [params.y0 zeros(size(params.y0,1), length(t)-p)];

if ~sys.tvp % simulate outputs with normal ARMAX parametrization
    % loop over all time steps
    for k=p:length(t)-1
        y(:, k+1) = sys.B_bar{1}*params.u(:, k+1);
        for i=1:p
            y(:, k+1) = y(:, k+1) + sys.A_bar{i,1} * y(:, k+1-i) ...
                + sys.B_bar{i+1} * params.u(:, k+1-i);
        end
    end
else % simulate outputs with time-varying parametrization
    y_init = reshape(params.y0, [],1);
    k=p;
    while k<=length(t)-1
        k_plus = k + p - 1;
        if ~isfield(sys, "A_tilde") || length(sys.A_tilde) < k+1 || isempty(sys.A_tilde(k+1))
            computeTVP(sys,k, 0:k_plus);
        end
        y_add = sys.A_tilde{k+1} * y_init;
        for i = 0:k_plus
            if k_plus -i <= length(t)-1
                y_add = y_add + sys.B_tilde{i+1,k+1} * params.u(:,k_plus + 1 -i);
            end
        end
        y(:, k+1:k_plus+1) = reshape(y_add,sys.nrOfOutputs,[]);
        k = k+p;
    end
    y = y(:, 1:length(t));
end
y=y';
end


% Auxiliary functions -----------------------------------------------------

function params = aux_uwv(sys,params,steps)

% check input, set default
if ~isfield(params,'u')
    params.u = zeros(sys.nrOfInputs,steps+1);
else
    % length always steps+1
    if size(params.u,2) == 1
        params.u = repmat(params.u,1,steps+1);
    elseif size(params.u,2) == steps && ~any(any(sys.B_bar{1}))
        params.u = [params.u, zeros(size(params.u,1),1)];
    elseif size(params.u,2) ~= steps+1
        throw(CORAerror('CORA:wrongFieldValue','params.u',...
            'has to be of size(2) of steps+1'));
    end
end
end

% ------------------------------ END OF CODE ------------------------------
