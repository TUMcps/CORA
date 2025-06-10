function [V,H,krylovOrder,Hlast] = priv_subspace_Krylov_individual_Jawecki(A,v,initKrylovOrder,options)
% priv_subspace_Krylov_individual_Jawecki - computes the Krylov subspace for a single
% vector given the accuracy to be achieved; the a-posteriori approach in
% equation 3.6 of [1]
%
%
% is used for tight error computation
%
% Syntax:
%    [V,H,krylovOrder,Hlast] = priv_subspace_Krylov_individual_Jawecki(A,v,initKrylovOrder,options)
%
% Inputs:
%    A - system matrix
%    nu_A - minimum eigenvalue of A
%    v - vector
%    initKrylovOrder - Krylov error that is first tested
%    options - reachability options
%              needs: krylovError, krylovStep and tFinal
%
% Outputs:
%    V - orthonormal basis
%    H - Hessenberg matrix
%    KrylovOrder - dimension of the reduced system
%    Hlast - last computed H_ij scalar of the Arnoldi iteration
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

% set precision for variable precison toolbox
precision = 34;

% compute norm of v
v_norm = norm(v);

% maximum allowed error
maxRelError = options.krylovError;

% initialize Krylov order and normalized error
krylovOrder = initKrylovOrder - options.krylovStep;
if krylovOrder <= 0
    krylovOrder = 1;
end

errorBound = inf;
dim = length(A);

% Krylov order should not be larger than dimension
if krylovOrder > dim
    krylovOrder = dim;
end

while (v_norm*errorBound > maxRelError) && (krylovOrder <= dim)
    
    % disp(errorBound);
    % increment Krylov order
    krylovOrder = krylovOrder + options.krylovStep;

    % perform Arnoldi iteration
    [V,H,Hlast,happyBreakdown] = arnoldi(A,v,krylovOrder);
    
    % sparsify
    V = sparse(V);
    H = sparse(H);

    % compute error if happy breakdown did not occur
    if ~happyBreakdown
        % compute/convert necessary parameters in mp
        tau = mp(Hlast,precision);
        
        gamma = 1;

        for i = 1:krylovOrder-1
            gamma = gamma*mp(H(i+1,i),precision);
        end

        % compute bound according to (Err_a) in [1]
        % for t we choose t_f since it leads to the largest error using this
        % timestep size
        % this line is computed seperately so the product isn't done in mp
        right_half = options.tFinal^krylovOrder/factorial(krylovOrder);
        errorBound = tau*gamma*right_half;

        if isnan(errorBound) % if error is not a number (NaN)
            errorBound = inf;
        end
    else
        % decrease Krylov order
        krylovOrder = length(V(1,:));
        break
    end
end

% ------------------------------ END OF CODE ------------------------------
