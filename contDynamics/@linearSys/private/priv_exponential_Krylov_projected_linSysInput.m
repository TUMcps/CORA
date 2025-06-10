function [R,R_Krylov] = priv_exponential_Krylov_projected_linSysInput(c_sys,g_sys,R0,dim_proj,options)
% priv_exponential_Krylov_projected_linSysInput - computes the overapproximation of the 
%    exponential of a system matrix up to a certain accuracy using a Krylov 
%    subspace; the subspace is already precomputed
%
% Syntax:
%    [R,R_Krylov] = priv_exponential_Krylov_projected_linSysInput(c_sys,g_sys,R0,C,options)
%
% Inputs:
%    c_sys - Linear system with sys.A = H_c, sys.B = V_c'
%    g_sys - Cell array of linear systems with sys{i}.A = H_g{i}, sys{i}.B = V_g{i}'
%    R0 - initial set
%    C - C matrix of linear system
%    options - reachability options
% 
%
% Outputs:
%    R - reachable set
%    R_Krylov - reachable set due to Krylov error
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Matthias Althoff, Maximilian Perschl
% Written:       02-November-2018
% Last update:   21-August-2024 (MP, replace stateflag with V/H as inputs)
%                25-April-2025 (MP, take linearSys as input to avoid recomputations)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Multiply previous reachable set with exponential matrix------------------
% obtain center and generators of previous reachable set
c = sparse(center(R0));
G = sparse(generators(R0));


c_norm = norm(c);
% check if center is zero
if isempty(c_sys)
    c_new = zeros(dim_proj,1);
else
    %Compute new center
    expMatrix = expm(c_sys.A*options.t);
    c_new = c_norm*c_sys.B'*expMatrix(:,1);
end


% preallocation
nrOfGens = size(G,2);
g_norm = zeros(nrOfGens,1);

% obtain generators using the Arnoldi iteration
G_new = zeros(dim_proj,size(G,2));

if nrOfGens > 0
    for iGen = 1:nrOfGens
    % parfor iGen = 1:nrOfGens
       g_norm(iGen) = norm(G(:,iGen));
        if g_norm(iGen) == 0
            G_new(:,iGen) = G(:,iGen);
        else
            %Compute new generator
            expMatrix = expm(g_sys{iGen}.A*options.t);
            G_new(:,iGen) = g_norm(iGen)*g_sys{iGen}.B'*expMatrix(:,1);
        end
        %iGen
    end
else
    G_new = []; %no generators 
end


if options.krylovError > 2*eps
    % Krylov error computation
    % +1 due to center
    error = options.krylovError * (size(R0.G,2) + 1);
    Krylov_interval = interval(-ones(length(c_new),1),ones(length(c_new),1))*error;
    R_Krylov = zonotope(Krylov_interval);

    % initial-state-solution zonotope
    R = zonotope([c_new,G_new,error*eye(length(c_new))]);
    %--------------------------------------------------------------------------
else
    R_Krylov = zonotope(0);
    R = zonotope([c_new,G_new]);
end

  
% ------------------------------ END OF CODE ------------------------------
