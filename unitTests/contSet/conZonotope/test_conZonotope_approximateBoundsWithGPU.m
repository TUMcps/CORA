function res = test_conZonotope_approximateBoundsWithGPU()
% test_conZonotope_approximateBoundsWithGPU - unit test function for 
%   GPU-based bound approximation of a batched constrained zonotopes.
%
% Syntax:
%    res = test_conZonotope_approximateBoundsWithGPU()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% References: 
%    [1] Koller, L. "Out of the Shadows: Exploring a Latent Space for 
%       Neural Network Verification". (2025) arXiv
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Lukas Koller
% Written:       19-November-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

rng('default')

% Specify the numerical tolerance.
tol = 1e-5;

% Specify number of dimensions, number of generators, and number of batch
% sizes.
n = 16; % Number of dimensions.
q = 32; % Number of generators.
bSz = 1; % Batch size.
p = 256; % Number of constraints.

% Sample a batch of random constraint zonotope.
cZs.c = rand([n bSz]);
cZs.G = rand([n q bSz]);
cZs.dr = rand([n bSz]);
cZs.A = rand([p q bSz]);
cZs.b = rand([p bSz]);

% Specify the options.
numUnionConst = 1; % Compute the intersection of all constraints.
options = nnHelper.validateNNoptions(struct());

% Approximate the bounds.
options.nn.conzonotope_bounding_method = 'fourier-motzkin';
options.nn.polytope_bound_approx_max_iter = 8;
[lFM,uFM,blFM,buFM] = conZonotope.approximateBoundsWithGPU(cZs,numUnionConst,options);

options.nn.conzonotope_bounding_method = 'dual-iter';
options.nn.conzonotope_bound_max_iter = 1000;
options.nn.conzonotope_bound_step_size = 1; 
[lD,uD,blD,buD] = conZonotope.approximateBoundsWithGPU(cZs,numUnionConst,options);

% Compute the exact bounds.
options.nn.conzonotope_bounding_method = 'exact';
[lE,uE,blE,buE] = conZonotope.approximateBoundsWithGPU(cZs,numUnionConst,options);

% Convert the inequality constraints to equality constraints.
cZeq = aux_2ConZonoWithEqConst(cZs);

% Check the computed bounds.
for i=1:bSz    
    % Construct the i-th constrained zonotope.
    cZi = conZonotope(cZeq.c(:,i),cZeq.G(:,:,i),cZeq.A(:,:,i),cZeq.b(:,i));

    % Compute the exact bounds using CORA.
    ivalcZi = interval(cZi);
    % Compare the bounds.
    assert(all(abs(lE(:,i) - ivalcZi.inf) <= tol)); % Check the lower bound.
    assert(all(abs(uE(:,i) - ivalcZi.sup) <= tol)); % Check the upper bound.
    
    % Check that the approximate bounds enclose the exact bounds.
    assert(all(lFM(:,i) <= lE(:,i) + tol)); % Check the lower bound.
    assert(all(uE(:,i) - tol <= uFM(:,i))); % Check the upper bound.

    % Check that the approximate bounds enclose the exact bounds.
    assert(all(lD(:,i) <= lE(:,i) + tol)); % Check the lower bound.
    assert(all(uE(:,i) - tol <= uD(:,i))); % Check the upper bound.

    % The 'dual-iter' bounds always have to be tighter.
    assert(all(lFM(:,i) <= lD(:,i) + tol)); % Check the lower bound.
    assert(all(uD(:,i) - tol <= uFM(:,i))); % Check the upper bound.

    % Check the bounds of the hypercube.
    assert(all(-1 <= blE(:,i) & blE(:,i) <= 1));
    assert(all(-1 <= buE(:,i) & buE(:,i) <= 1));
    assert(all(-1 <= blFM(:,i) & blFM(:,i) <= 1));
    assert(all(-1 <= buFM(:,i) & buFM(:,i) <= 1));

    % Construct the constrained hypercube.
    Ai = double(gather(cZs.A(:,:,i)));
    bi = double(gather(cZs.b(:,i)));
    cPi = interval(-ones([q 1]),ones([q 1])) & polytope(Ai,bi);
    % Compute the exact bounds of the hypercube.
    ivalcPi = interval(cPi);
    % Compare the exact bounds.
    assert(all(withinTol(blE(:,i),ivalcPi.inf,tol))); % Check the lower bound.
    assert(all(withinTol(buE(:,i),ivalcPi.sup,tol))); % Check the upper bound.
    % Check that the approximated bounds enclose the exact bounds.
    assert(all(blFM(:,i) <= ivalcPi.inf + tol)); % Check the lower bound.
    assert(all(ivalcPi.sup - tol <= buFM(:,i))); % Check the upper bound.
end

% Set the result.
res = true;

end


% Auxiliary functions -----------------------------------------------------

function cZeq = aux_2ConZonoWithEqConst(cZineq)
    % Extract parameters of the constraint zonotope.
    c = double(gather(cZineq.c));
    G = double(gather(cZineq.G));
    dr = double(gather(cZineq.dr));
    A = double(gather(cZineq.A));
    b = double(gather(cZineq.b));

    % We convert the inequality constraints to equality constraints by 
    % adding a slack variable.

    % Obtain number of dimensions, generators, and batch size.
    [n,~,bSz] = size(G);
    % Obtain number of constraints.
    [p,~] = size(A);

    cZeq.c = c;
    % Add the radius to the generators.
    if any(dr ~= 0,'all')
        G = cat(2,G,permute(dr,[1 3 2]).*eye(n));
        A = cat(2,A,zeros([p n bSz]));
    end
    % Add a slack variable.
    cZeq.G = cat(2,G,zeros([n p bSz]));
    % Compute scale for the slack variable.
    s = 1/2*(sum(abs(A),2) + permute(b,[1 3 2]));
    cZeq.A = cat(2,A,eye(p).*s);
    % Compensate for the slack variable.
    cZeq.b = b - s(:,:);
end

% ------------------------------ END OF CODE ------------------------------
