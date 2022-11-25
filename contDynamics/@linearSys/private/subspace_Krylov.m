function [V_c,H_c,V_g,H_g] = subspace_Krylov(obj,Z,options)
% subspace_Krylov - computes the Krylov subspaces for the center and
% generators of a zonotpe
%
% Syntax:  
%    [obj] = exponential_Krylov(obj,options)
%
% Inputs:
%    obj - linearSys object
%    Z - zonotope
%    options - reachability options
%
% Outputs:
%    V_c - orthonormal basis of center
%    H_c - Hessenberg matrix of center
%    V_g - orthonormal basis of generators
%    H_g - Hessenberg matrix of generators
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
c = sparse(center(Z));
G = sparse(generators(Z));

% init Krylov order
KrylovOrder = 1;

% change obj.A to fit paper
A = obj.A;

% minimum eigenvalue
%nu_A_test = eigs((-A-A')/2, 1, 'sa'); % lambda_min(A+A*/2) %<-- for Wang approx.
nu_A = eigs((A+A')/2, 1, 'lm'); % lambda_max(A+A*/2) %<-- for Jia approx.

% check if center is zero
c_norm = norm(c);
if c_norm == 0
    V_c = [];
    H_c = [];
else
    %Arnoldi; Krylov order is passed on to computation of next generator
    [V_c,H_c,KrylovOrder] = subspace_Krylov_individual(A,nu_A,c,KrylovOrder,options);
    % change results back \dot{x} = Ax (Krylov papers use \dot{x} = -Ax);
    % V remains unchanged
    %H_c = -H_c;
end
    
% preallocation
nrOfGens = length(G(1,:));
V_g = cell(nrOfGens,1);
H_g = cell(nrOfGens,1);

% obtain generators using the Arnoldi iteration
if nrOfGens>0
    for iGen = 1:nrOfGens
        g_norm = norm(G(:,iGen));
        if g_norm == 0
            V_g{iGen} = [];
            H_g{iGen} = [];
        else
            %Arnoldi
            [V_g{iGen},H_g{iGen},KrylovOrder] = subspace_Krylov_individual(A,nu_A,G(:,iGen),KrylovOrder,options);
            % change results back \dot{x} = Ax (Krylov papers use \dot{x} = -Ax);
            % V remains unchanged
            %H_g{iGen} = -H_g{iGen};
        end
        iGen
        KrylovOrder
    end
else
    V_g = []; %no generators 
    H_g = []; %no generators 
end

%------------- END OF CODE --------------