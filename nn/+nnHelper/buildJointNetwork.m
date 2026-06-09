function [nn, X0_reduced, netInfo] = buildJointNetwork(nn_f, nn_g, X0_joint, info)
% buildJointNetwork - assemble a joint network [f(x_f); g(x_g)] for
%    multi-network VNN-COMP benchmarks (monotonic/isomorphic ACAS Xu): detect
%    the input couplings X_f[i]==X_g[i], build selector matrices Sf/Sg and the
%    concat composite network, and project the input set onto the free vars.
%
% Syntax:
%    [nn, X0_reduced, netInfo] = nnHelper.buildJointNetwork(nn_f, nn_g, X0_joint, info)
%
% Inputs:
%    nn_f - first sub-network (neuralNetwork)
%    nn_g - second sub-network (neuralNetwork)
%    X0_joint - cell {polytope} on the joint input space (dim = nf+ng)
%    info - parser info struct (.inputDims, .outputDims, .networks)
%
% Outputs:
%    nn - joint neuralNetwork (single concat nnCompositeLayer)
%    X0_reduced - cell {interval|polytope} on the free-variable input space
%    netInfo - metadata for counterexample formatting (.netNames, .Sf, .Sg,
%              .E, .coupledF, .coupledG, .uncoupledG, .nFree, .nfOut, .ngOut)
%
% Other m-files required: none
% Subfunctions: aux_findEqualityCouplings
% MAT-files required: none
%
% See also: prepare_instance, neuralNetwork, nnCompositeLayer, vnnlib2cora

% Authors:       Benedikt Kellner
% Written:       03-June-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% X0_joint is a cell array {polytope} from vnnlib2cora; extract the polytope.
if iscell(X0_joint)
    X0_poly = X0_joint{1};
else
    X0_poly = X0_joint;
end

% only two sub-networks (joint [f; g]) are supported
if numel(info.inputDims) ~= 2
    throw(CORAerror('CORA:notSupported', ...
        sprintf(['buildJointNetwork supports exactly 2 sub-networks, ' ...
        'got %d.'], numel(info.inputDims))));
end

nf = info.inputDims(1);
ng = info.inputDims(2);

% Find which g-input dims are equality-coupled to f-input dims.
[coupledF, coupledG] = aux_findEqualityCouplings(X0_poly, nf, ng);

% Free variables: all of X_f (nf dims) + uncoupled dims of X_g.
uncoupledG = setdiff(1:ng, coupledG);
nFree = nf + numel(uncoupledG);

% Selector matrices  x_f = Sf*v,  x_g = Sg*v  (each is n*×nFree).
Sf = [eye(nf), zeros(nf, numel(uncoupledG))];

Sg = zeros(ng, nFree);
for k = 1:numel(coupledG)
    Sg(coupledG(k), coupledF(k)) = 1;     % X_g[coupled] = v[coupledF partner]
end
for k = 1:numel(uncoupledG)
    Sg(uncoupledG(k), nf + k) = 1;        % X_g[uncoupled] = v[nf+k]
end

% Build each computation path: affine input selector + network layers.
path_f = [{nnLinearLayer(Sf, zeros(nf,1))}, nn_f.layers(:)'];
path_g = [{nnLinearLayer(Sg, zeros(ng,1))}, nn_g.layers(:)'];
nn = neuralNetwork({nnCompositeLayer({path_f; path_g}, 'concat')});
% set inputSize explicitly (not inferable when layer 1 is a composite layer)
nn.setInputSize([nFree, 1]);

% project the joint polytope onto the free variables: x_joint = E * v
n = nf + ng;
E = zeros(n, nFree);
E(1:nf, 1:nf) = eye(nf);
for k = 1:numel(coupledG)
    E(nf + coupledG(k), coupledF(k)) = 1;
end
for k = 1:numel(uncoupledG)
    E(nf + uncoupledG(k), nf + k) = 1;
end

A_proj = X0_poly.A * E;
b_proj = X0_poly.b;

% Remove trivial rows (sum of |coefficients| == 0, constraint 0 <= b >= 0).
tol = 1e-8;
rowNnz = sum(abs(A_proj) > tol, 2);
keep   = rowNnz > 0;
A_nt   = A_proj(keep, :);
b_nt   = b_proj(keep);

% Determine if the remaining constraints form a box.
if all(rowNnz(keep) == 1)
    % All constraints are axis-aligned: extract interval bounds directly.
    lb = -inf(nFree, 1);
    ub =  inf(nFree, 1);
    for r = 1:size(A_nt, 1)
        col = find(abs(A_nt(r,:)) > tol);
        v   = A_nt(r, col);
        if v > 0
            ub(col) = min(ub(col), b_nt(r) / v);
        else
            lb(col) = max(lb(col), b_nt(r) / v);
        end
    end
    X0_reduced = {interval(lb, ub)};
else
    % Non-box (e.g., monotonic benchmark has residual X_g0 <= X_f0 coupling).
    X0_reduced = {polytope(A_nt, b_nt)};
end

% Build metadata struct for multi-net counterexample formatting / inspection.
netInfo.netNames = {info.networks(1).name, info.networks(2).name};
netInfo.Sf = Sf;
netInfo.Sg = Sg;
netInfo.E = E;
netInfo.coupledF = coupledF;
netInfo.coupledG = coupledG;
netInfo.uncoupledG = uncoupledG;
netInfo.nFree = nFree;
netInfo.nfOut = info.outputDims(1);
netInfo.ngOut = info.outputDims(2);

end


% Auxiliary functions -----------------------------------------------------

function [coupledF, coupledG] = aux_findEqualityCouplings(X0, nf, ng)
% paired indices with X_f[coupledF(k)]==X_g[coupledG(k)], from the polytope
% equality rows (.Ae/.be). X_f is columns 1..nf, X_g columns nf+1..nf+ng.

Ae = X0.Ae;
be = X0.be;
m  = size(Ae, 1);

coupledF = zeros(1, m);
coupledG = zeros(1, m);
tol = 1e-9;

for j = 1:m
    nz = find(abs(Ae(j,:)) > tol);
    fi = nz(nz <= nf);          % nonzero column(s) in the X_f block
    gi = nz(nz >  nf) - nf;     % nonzero column(s) in the X_g block
    % a coupling pairs one X_f dim with one X_g dim, with be == 0
    if numel(nz) == 2 && isscalar(fi) && isscalar(gi) && abs(be(j)) <= tol
        coupledF(j) = fi;
        coupledG(j) = gi;
    else
        % reject any other equality (same-network, scaled, or nonzero rhs)
        throw(CORAerror('CORA:notSupported', sprintf( ...
            ['buildJointNetwork: input equality constraint %d is not a ' ...
             'simple X_f[i]==X_g[j] coupling.'], j)));
    end
end

end

% ------------------------------ END OF CODE ------------------------------
