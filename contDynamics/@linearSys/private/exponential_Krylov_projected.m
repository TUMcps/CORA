function [R,R_Krylov] = exponential_Krylov_projected(linsys,R0,options,stateFlag)
% exponential_Krylov_projected - computes the overapproximation of the 
%    exponential of a system matrix up to a certain accuracy using a Krylov 
%    subspace; the subspace is already precomputed
%
% Syntax:
%    [R,R_Krylov] = exponential_Krylov_projected(linsys,R0,options,stateFlag)
%
% Inputs:
%    linsys - linearSys object
%    R0 - initial set
%    options - reachability options
%    stateFlag - true if computation is for the state solution
%
% Outputs:
%    R - reachable set
%    R_Krylov - rechable set due to Krylov error
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       02-November-2018
% Last update:   02-June-2020
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Multiply previous reachable set with exponential matrix

% obtain center and generators of previous reachable set
c = sparse(R0.c);
G = sparse(R0.G);

% retrieve V_c and H_c from Arnoldi
if stateFlag % for state solution
    V_c_proj = linsys.krylov.state.V_c_proj;
    V_g_proj = linsys.krylov.state.V_g_proj;
    H_c = linsys.krylov.state.H_c;
    H_g = linsys.krylov.state.H_g;
    options.t = options.t + options.timeStep; % state solution starts from t=0; input solution from t\in[0,timeStep]
else % for input solution
    V_c_proj = linsys.krylov.input.V_c_proj;
    V_g_proj = linsys.krylov.input.V_g_proj;
    H_c = linsys.krylov.input.H_c;
    H_g = linsys.krylov.input.H_g;
end

% check if center is zero
c_norm = norm(c);
if c_norm == 0
    if ~isscalar(linsys.C) || linsys.C~=1
        dim_proj = length(linsys.C(:,1));
    else
        dim_proj = length(c);
    end
    c_new = zeros(dim_proj,1);
else
    %Compute new center
    expMatrix = expm(H_c*options.t);
    c_new = c_norm*V_c_proj*expMatrix(:,1);
end


% preallocation
nrOfGens = length(G(1,:));
g_norm = zeros(nrOfGens,1);

% obtain generators using the Arnoldi iteration
if nrOfGens>0
    parfor iGen = 1:nrOfGens
    %for iGen = 1:nrOfGens
        g_norm(iGen) = norm(G(:,iGen));
        if g_norm(iGen) == 0
            G_new(:,iGen) = G(:,iGen);
        else
            %Compute new generator
            expMatrix = expm(H_g{iGen}*options.t);
            G_new(:,iGen) = g_norm(iGen)*V_g_proj{iGen}*expMatrix(:,1);
        end
        %iGen
    end
else
    G_new = []; %no generators 
end

if options.krylovError > 2*eps
    % Krylov error computation
    error_normalized = options.krylovError;
    err = error_normalized*options.t;
    Krylov_interval = interval(-1,1)*ones(length(c_new),1)*err;
    R_Krylov = zonotope(Krylov_interval);

    % initial-state-solution zonotope
    R = zonotope(c_new,G_new,err*eye(length(c_new)));
else
    R_Krylov = 0;
    R = zonotope(c_new,G_new);
end
  
% ------------------------------ END OF CODE ------------------------------
