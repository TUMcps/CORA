function [V_c,H_c,V_g,H_g,Hlast_c] = priv_subspace_Krylov_jaweckiBound(A,Z,options)
% priv_subspace_Krylov_jaweckiBound - computes the Krylov subspaces for the center and
% generators of a zonotpe
% however, the error bound is computed as given in [1]
%
% Syntax:
%    [V_c,H_c,V_g,H_g,Hlast_c] = priv_subspace_Krylov_jaweckiBound(A,Z,options)
%
% Inputs:
%    A - state matrix of linear system
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
% References:
%    [1] Computable upper error bounds for Krylov approximations to
%    matrix exponentials and associated phi-functions, Jawecki et al, 2019
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Maximilian Perschl
% Written:       25-April-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Multiply previous reachable set with exponential matrix------------------
% obtain center and generators of previous reachable set
c = sparse(center(Z));
G = sparse(generators(Z));

% init Krylov order
KrylovOrder = 15;

% check if center is zero
c_norm = norm(c);
if c_norm == 0
    V_c = [];
    H_c = [];
else
    %Arnoldi; Krylov order is passed on to computation of next generator
    [V_c,H_c,KrylovOrder,Hlast_c] = priv_subspace_Krylov_individual_Jawecki(A,c,KrylovOrder,options);
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
            [V_g{iGen},H_g{iGen},KrylovOrder] = priv_subspace_Krylov_individual_Jawecki(A,G(:,iGen),KrylovOrder,options);
        end
%         debug output
%         iGen
        % KrylovOrder
    end
else
    V_g = []; %no generators 
    H_g = []; %no generators 
end

end

% ------------------------------ END OF CODE ------------------------------
