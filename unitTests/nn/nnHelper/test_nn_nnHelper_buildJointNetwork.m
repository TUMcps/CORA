function res = test_nn_nnHelper_buildJointNetwork()
% test_nn_nnHelper_buildJointNetwork - test the joint-network builder:
%    coupling detection, selector matrices, input-set projection
%    (box -> interval, non-box -> polytope), and evaluation equivalence
%    nn.evaluate(v) == [f(Sf*v); g(Sg*v)] for points and intervals
%
% Syntax:
%    res = test_nn_nnHelper_buildJointNetwork()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nnHelper.buildJointNetwork, prepare_instance

% Authors:       Benedikt Kellner
% Written:       03-June-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Reproducible synthetic networks.
rng(1234);

% =====================================================================
% Case A: fully coupled  X_f[i] == X_g[i]  ->  interval X0
% =====================================================================
nf = 3; ng = 3; mf = 2; mg = 4;
nn_f = aux_makeNet(nf, mf);
nn_g = aux_makeNet(ng, mg);
info = aux_makeInfo(nf, ng, mf, mg);

lo = [-1; -2; 0]; hi = [1; 0.5; 3];
[Cb, db] = aux_boxRows(lo, hi, nf, ng, 'f');
[Ce, de] = aux_equalityRows([1 2 3], [1 2 3], nf, ng);
X0 = {polytope(Cb, db, Ce, de)};

[nn, X0red, netInfo] = nnHelper.buildJointNetwork(nn_f, nn_g, X0, info);

% Coupling detection.
assert(isequal(sort(netInfo.coupledF), [1 2 3]), 'A: coupledF wrong');
assert(isequal(sort(netInfo.coupledG), [1 2 3]), 'A: coupledG wrong');
assert(isempty(netInfo.uncoupledG), 'A: expected no uncoupled g-dims');
assert(netInfo.nFree == 3, 'A: nFree should be 3');

% Selector matrices: fully coupled -> Sf = Sg = eye(3).
assert(isequal(netInfo.Sf, eye(3)), 'A: Sf should be eye(3)');
assert(isequal(netInfo.Sg, eye(3)), 'A: Sg should be eye(3)');

% Box couplings -> interval input set with the original X_f bounds.
assert(isscalar(X0red) && isa(X0red{1}, 'interval'), 'A: X0 not interval');
assert(dim(X0red{1}) == 3, 'A: X0 dim should be 3');
assert(all(withinTol(infimum(X0red{1}), lo, 1e-9)), 'A: X0 inf wrong');
assert(all(withinTol(supremum(X0red{1}), hi, 1e-9)), 'A: X0 sup wrong');

% Network metadata.
assert(isequal(netInfo.netNames, {'f','g'}), 'A: netNames wrong');
assert(netInfo.nfOut == mf && netInfo.ngOut == mg, 'A: out dims wrong');
assert(isa(nn, 'neuralNetwork'), 'A: nn not a neuralNetwork');

% Evaluation equivalence.
aux_assertPointEquiv(nn, nn_f, nn_g, netInfo, 'A');
aux_assertIntervalEquiv(nn, nn_f, nn_g, netInfo, 'A');

% =====================================================================
% Case B: partial coupling + uncoupled inequality  ->  polytope X0
%    couple dims 2,3 (X_f[i]==X_g[i]); dim 1 free with X_g[1] <= X_f[1]
% =====================================================================
nf = 3; ng = 3; mf = 2; mg = 2;
nn_f = aux_makeNet(nf, mf);
nn_g = aux_makeNet(ng, mg);
info = aux_makeInfo(nf, ng, mf, mg);

[Cb, db]   = aux_boxRows([-1;-1;-1], [1;1;1], nf, ng, 'f'); % box on X_f
[Cbg, dbg] = aux_boxRows(-2, 2, nf, ng, 'g');               % box on X_g[1]
[Ce, de]   = aux_equalityRows([2 3], [2 3], nf, ng);        % X_f[i]==X_g[i]
% Inequality X_g[1] <= X_f[1]  ->  -X_f[1] + X_g[1] <= 0  (NOT an equality)
n = nf + ng;
rIneq = zeros(1, n); rIneq(1) = -1; rIneq(nf+1) = 1;
X0 = {polytope([Cb; Cbg; rIneq], [db; dbg; 0], Ce, de)};

[nn, X0red, netInfo] = nnHelper.buildJointNetwork(nn_f, nn_g, X0, info);

assert(isequal(sort(netInfo.coupledF), [2 3]), 'B: coupledF wrong');
assert(isequal(sort(netInfo.coupledG), [2 3]), 'B: coupledG wrong');
assert(isequal(netInfo.uncoupledG, 1), 'B: uncoupledG should be {1}');
assert(netInfo.nFree == 4, 'B: nFree should be 4');
% The X_g[1] <= X_f[1] inequality must NOT be mistaken for a coupling.
assert(~ismember(1, netInfo.coupledG), 'B: inequality misread as coupling');

% Non-box residual -> polytope of the free-variable dimension.
assert(isscalar(X0red) && isa(X0red{1}, 'polytope'), 'B: X0 not polytope');
assert(dim(X0red{1}) == 4, 'B: X0 dim should be 4');

% Selector shapes: Sf in R^{3x4}, Sg in R^{3x4}.
assert(isequal(size(netInfo.Sf), [3 4]), 'B: Sf shape wrong');
assert(isequal(size(netInfo.Sg), [3 4]), 'B: Sg shape wrong');

aux_assertPointEquiv(nn, nn_f, nn_g, netInfo, 'B');
aux_assertIntervalEquiv(nn, nn_f, nn_g, netInfo, 'B');

% =====================================================================
% Case C: no couplings (two independent nets)  ->  interval X0
% =====================================================================
nf = 2; ng = 2; mf = 3; mg = 1;
nn_f = aux_makeNet(nf, mf);
nn_g = aux_makeNet(ng, mg);
info = aux_makeInfo(nf, ng, mf, mg);

[Cf, df] = aux_boxRows([-1;-1], [1;1], nf, ng, 'f');
[Cg, dg] = aux_boxRows([-3;0],  [3;2], nf, ng, 'g');
X0 = {polytope([Cf; Cg], [df; dg])};

[nn, X0red, netInfo] = nnHelper.buildJointNetwork(nn_f, nn_g, X0, info);

assert(isempty(netInfo.coupledF) && isempty(netInfo.coupledG), 'C: spurious coupling');
assert(isequal(sort(netInfo.uncoupledG), [1 2]), 'C: uncoupledG wrong');
assert(netInfo.nFree == 4, 'C: nFree should be 4');
assert(isequal(netInfo.Sf, [eye(2), zeros(2,2)]), 'C: Sf wrong');
assert(isequal(netInfo.Sg, [zeros(2,2), eye(2)]), 'C: Sg wrong');
assert(isa(X0red{1}, 'interval') && dim(X0red{1}) == 4, 'C: X0 not 4-D interval');

aux_assertPointEquiv(nn, nn_f, nn_g, netInfo, 'C');
aux_assertIntervalEquiv(nn, nn_f, nn_g, netInfo, 'C');

% test completed
res = true;

end


% Auxiliary functions -----------------------------------------------------

function nn = aux_makeNet(nin, nout)
% Small 2-layer ReLU network  nin -> 5 -> nout.
W1 = rand(5, nin) - 0.5; b1 = rand(5, 1) - 0.5;
W2 = rand(nout, 5) - 0.5; b2 = rand(nout, 1) - 0.5;
nn = neuralNetwork({nnLinearLayer(W1, b1), nnReLULayer(), nnLinearLayer(W2, b2)});
nn.setInputSize([nin 1]);
end

function info = aux_makeInfo(nf, ng, mf, mg)
info.inputDims  = [nf ng];
info.outputDims = [mf mg];
info.networks(1).name = 'f';
info.networks(2).name = 'g';
end

function [C, d] = aux_boxRows(lo, hi, nf, ng, which)
% Axis-aligned box  lo <= X_<which> <= hi  on the joint (nf+ng)-D space.
n = nf + ng;
if strcmp(which, 'f'), off = 0; else, off = nf; end
k = numel(lo);
C = []; d = [];
for i = 1:k
    r = zeros(1, n); r(off+i) =  1; C = [C; r]; d = [d;  hi(i)]; %#ok<AGROW>
    r = zeros(1, n); r(off+i) = -1; C = [C; r]; d = [d; -lo(i)]; %#ok<AGROW>
end
end

function [Ce, de] = aux_equalityRows(fi, gi, nf, ng)
% Equality couplings  X_f[fi(k)] == X_g[gi(k)]  as native equality rows
% (Ae*x == be), matching how the VNN-LIB 2.0 parser emits them into the
% joint input polytope's equality constraints.
n  = nf + ng;
Ce = zeros(numel(fi), n);
de = zeros(numel(fi), 1);
for k = 1:numel(fi)
    Ce(k, fi(k)) = 1; Ce(k, nf + gi(k)) = -1;   % X_f[fi] - X_g[gi] == 0
end
end

function aux_assertPointEquiv(nn, nn_f, nn_g, netInfo, tag)
% nn.evaluate(v) == [f(Sf*v); g(Sg*v)] for a batch of random points.
V = 2 * rand(netInfo.nFree, 64) - 1;
Y = nn.evaluate(V);
Yref = [nn_f.evaluate(netInfo.Sf * V); nn_g.evaluate(netInfo.Sg * V)];
assert(all(abs(Y - Yref) < 1e-8, 'all'), ...
    sprintf('%s: point evaluation not equivalent to [f(Sf*v); g(Sg*v)]', tag));
end

function aux_assertIntervalEquiv(nn, nn_f, nn_g, netInfo, tag)
% nn.evaluate(Iv) == [f(Sf*Iv); g(Sg*Iv)] for an input interval.
Iv = interval(-0.3 * ones(netInfo.nFree, 1), 0.4 * ones(netInfo.nFree, 1));
Iy = nn.evaluate(Iv);
Iyref = [nn_f.evaluate(netInfo.Sf * Iv); nn_g.evaluate(netInfo.Sg * Iv)];
assert(all(withinTol(infimum(Iy), infimum(Iyref), 1e-7)), ...
    sprintf('%s: interval lower bound mismatch', tag));
assert(all(withinTol(supremum(Iy), supremum(Iyref), 1e-7)), ...
    sprintf('%s: interval upper bound mismatch', tag));
end

% ------------------------------ END OF CODE ------------------------------
