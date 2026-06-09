function res = testnn_vnnlib2cora_multinet()
% testnn_vnnlib2cora_multinet - check that the VNN-LIB 2.0 parser emits
%    cross-variable couplings (X_f[i]==X_g[j]) as polytope equality rows
%    (.Ae/.be), while box bounds, value-fixing equalities and residual
%    inequalities stay as inequality rows (.A/.b)
%
% Syntax:
%    res = testnn_vnnlib2cora_multinet()
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
% See also: vnnlib2cora, nnHelper.buildJointNetwork

% Authors:       Benedikt Kellner
% Written:       05-June-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% A minimal two-network spec: inputs X_f[0..1], X_g[0..1].
%   - X_f[0] == X_g[0]            -> coupling (equality)
%   - X_f[1] >= X_g[1]            -> residual inequality (NOT a coupling)
%   - X_f[1] == 0.5              -> value-fixing equality (NOT a coupling)
%   - box bounds on the rest      -> inequalities
spec = [ ...
    "(vnnlib-version <2.0>)" newline ...
    "(declare-network f (declare-input X_f real [2]) (declare-output Y_f real [2]))" newline ...
    "(declare-network g (equal-to f) (declare-input X_g real [2]) (declare-output Y_g real [2]))" newline ...
    "(assert (and (<= X_f[0] 1.0) (>= X_f[0] -1.0)))" newline ...
    "(assert (== X_f[1] 0.5))" newline ...
    "(assert (== X_f[0] X_g[0]))" newline ...
    "(assert (and (>= X_f[1] X_g[1]) (>= X_g[1] -1.0)))" newline ...
    "(assert (<= Y_f[0] Y_g[0]))" newline ];

specFile = [tempname '.vnnlib'];
cleanup = onCleanup(@() aux_safeDelete(specFile));
fid = fopen(specFile, 'w');
assert(fid > 0, 'could not open temp vnnlib file');
fwrite(fid, char(strjoin(spec, '')));
fclose(fid);

[X0, ~, info] = vnnlib2cora(specFile);

% --- parser metadata ---
assert(numel(info.networks) == 2, 'expected two declared networks.');
assert(isequal(info.inputDims, [2 2]), 'inputDims should be [2 2].');

% --- the joint input set is a polytope (non-box, coupled) ---
assert(isscalar(X0), 'expected a single input branch.');
P = X0{1};
assert(isa(P, 'polytope'), 'multi-network input set should be a polytope.');

% --- exactly ONE native equality row: the X_f[0]==X_g[0] coupling ---
assert(size(P.Ae, 1) == 1, ...
    sprintf('expected exactly 1 equality row, got %d.', size(P.Ae,1)));
% totalIn = 4 columns: [X_f0 X_f1 X_g0 X_g1]; coupling pairs col 1 and col 3.
nz = find(abs(P.Ae(1,:)) > 1e-9);
assert(isequal(nz(:)', [1 3]), 'coupling should pair X_f[0] (col1) and X_g[0] (col3).');
assert(abs(P.Ae(1,1)) == 1 && abs(P.Ae(1,3)) == 1, 'coupling coefficients should be +/-1.');
assert(sign(P.Ae(1,1)) == -sign(P.Ae(1,3)), 'coupling must have opposite signs.');
assert(abs(P.be(1)) <= 1e-9, 'coupling rhs should be 0.');

% --- the residual inequality and value-fixing equality stay in .A ---
assert(any(abs(P.A(:,2)) > 1e-9), 'value-fixing equality X_f[1]==0.5 should be in .A.');
% residual X_f[1] >= X_g[1] couples cols 2 and 4 in an inequality row of A.
hasResidual = any(sum(abs(P.A(:, [2 4])) > 1e-9, 2) == 2 & ...
                  sum(abs(P.A) > 1e-9, 2) == 2);
assert(hasResidual, 'residual X_f[1] >= X_g[1] inequality should be in .A.');

% --- and the interval enclosure still accounts for the equality ---
I = interval(P);
assert(dim(I) == 4, 'interval enclosure should be 4-D.');

% test completed
res = true;

end


% Auxiliary functions -----------------------------------------------------

function aux_safeDelete(f)
if isfile(f), delete(f); end
end

% ------------------------------ END OF CODE ------------------------------
