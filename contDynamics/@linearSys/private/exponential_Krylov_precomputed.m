function [R,R_Krylov] = exponential_Krylov_precomputed(obj,R0,options,stateFlag)
% exponential_krylov_precomputed - computes the overapproximation of the 
% exponential of a system matrix up to a certain accuracy using a Krylov 
% subspace; the subspace is already precomputed
%
% Syntax:  
%    [R,R_Krylov] = exponential_Krylov_precomputed(obj,options)
%
% Inputs:
%    obj - linearSys object
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
% See also: ---

% Author:       Matthias Althoff
% Written:      02-November-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% Multiply previous reachable set with exponential matrix------------------
% obtain center and generators of previous reachable set
c = sparse(center(R0));
G = sparse(generators(R0));

% retrieve V_c and H_c from Arnoldi
if stateFlag % for state solution
    V_c = obj.krylov.state.V_c;
    H_c = obj.krylov.state.H_c;
    V_g = obj.krylov.state.V_g;
    H_g = obj.krylov.state.H_g;
else % for input solution
    V_c = obj.krylov.state.V_c;
    H_c = obj.krylov.state.H_c;
    V_g = obj.krylov.state.V_g;
    H_g = obj.krylov.state.H_g;
end

% check if center is zero
c_norm = norm(c);
if c_norm == 0
    c_new = c;
else
    %Compute new center
    expMatrix = expm(H_c*options.t);
    c_new = c_norm*V_c*expMatrix(:,1);
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
            G_new(:,iGen) = g_norm(iGen)*V_g{iGen}*expMatrix(:,1);
        end
        %iGen
    end
else
    G_new = []; %no generators 
end

if options.krylovError > 2*eps
    % Krylov error computation
    error_normalized = obj.krylov.errorBound_normalized;
    error = norm(Rinit)*error_normalized*options.t;
    Krylov_interval = interval(-1,1)*ones(dimension(obj),1)*error;
    R_Krylov = zonotope(Krylov_interval);

    % initial-state-solution zonotope
    R = zonotope([c_new,G_new,error*eye(dimension(obj))]);
    %--------------------------------------------------------------------------
else
    R_Krylov = 0;
    R = zonotope([c_new,G_new]);
end

  

%------------- END OF CODE --------------