function [R,R_Krylov] = exponential_Krylov(obj,R,options)
% exponential_krylov - computes the overapproximation of the exponential of a system 
% matrix up to a certain accuracy using a Krylov subspace
%
% Syntax:  
%    [obj] = exponential_Krylov(obj,options)
%
% Inputs:
%    obj - linearSys object
%    linRed - reduced linear system
%    options - reachability options
%
% Outputs:
%    obj - linearSys object
%    linRed - reduced linear system
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
% Written:      03-March-2017 
% Last update:  25-October-2018
% Last revision:---

%------------- BEGIN CODE --------------

% Multiply previous reachable set with exponential matrix------------------
% obtain center and generators of previous reachable set
c = sparse(center(R));
G = sparse(generators(R));

% check if center is zero
c_norm = norm(c);
if c_norm == 0
    c_new = c;
else
    %Arnoldi
    [V_c,H_c] = arnoldi(obj,c,options.redDim);
    
    %Compute new center
    expMatrix = expm(H_c*options.timeStep);
    c_new = c_norm*V_c*expMatrix(:,1);
    
end


% preallocation
nrOfGens = length(G(1,:));
g_norm = zeros(nrOfGens,1);
V_g = cell(nrOfGens,1);
H_g = cell(nrOfGens,1);

% obtain generators using the Arnoldi iteration
if nrOfGens>0
    parfor iGen = 1:nrOfGens
    %for iGen = 1:nrOfGens
        g_norm(iGen) = norm(G(:,iGen));
        if g_norm(iGen) == 0
            G_new(:,iGen) = G(:,iGen);
        else
            %Arnoldi
            [V_g{iGen},H_g{iGen}] = arnoldi(obj,G(:,iGen),options.redDim);

            %Compute new generator
            expMatrix = expm(H_g{iGen}*options.timeStep);
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
    error = norm(Rinit)*error_normalized;
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