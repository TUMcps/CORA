function res = test_nn_converter_vnnlib2cora_v2_multinet()
% test_nn_converter_vnnlib2cora_v2_multinet - multi-network VNN-LIB 2.0
%    features: (declare-network ...) blocks, (isomorphic-to ...), joint
%    input/output spaces, cross-network linear constraints.
%
% Syntax:
%    res = test_nn_converter_vnnlib2cora_v2_multinet()
%
% Inputs:
%    -
%
% Outputs:
%    res - true on success
%
% See also: vnnlib2cora, vnnlib2cora_v2

% Authors:       Benedikt Kellner
% Written:       19-April-2026
% Last update:   08-May-2026 (BK, isomorphism type field assertions; chained equivalence error)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

file = [CORAROOT '/models/Cora/nn/unitTests/vnnlib/v2/acasxu_equivalence.vnnlib'];
[X0, spec, info] = vnnlib2cora(file);

% Header
assert(strcmp(info.version, '2.0'));

% Two networks declared, second is isomorphic to the first
assert(numel(info.networks) == 2);
assert(strcmp(info.networks(1).name, 'f'));
assert(strcmp(info.networks(2).name, 'g'));
assert(isempty(info.networks(1).isomorphicTo));
assert(strcmp(info.networks(2).isomorphicTo, 'f'));

assert(numel(info.isomorphism) == 1);
assert(strcmp(info.isomorphism(1).source, 'g'));
assert(strcmp(info.isomorphism(1).target, 'f'));
assert(strcmp(info.isomorphism(1).type, 'isomorphic'));
assert(isempty(info.networks(2).equalTo));

% Joint input/output dimensions: 5 + 5
assert(isequal(info.inputDims,  [5 5]));
assert(isequal(info.outputDims, [5 5]));
assert(info.totalIn == 10 && info.totalOut == 10);

% X0: a polytope on the 10-D joint input space with the equality coupling
% X_f[i] == X_g[i] baked in (so the polytope has both <= and >= for each
% pair, which means it cannot be downcast to an interval -> stays polytope)
assert(isscalar(X0));
assert(isa(X0{1}, 'polytope'));
assert(dim(X0{1}) == 10);

% spec: scalar specification on the 10-D joint output space
assert(isscalar(spec));
assert(dim(spec.set) == 10);

% --- infix comparison form: (lhs op rhs) instead of prefix (op lhs rhs) ---
file_infix = [CORAROOT '/models/Cora/nn/unitTests/vnnlib/v2/infix_comparison.vnnlib'];
[X0_infix, spec_infix, info_infix] = vnnlib2cora(file_infix);
assert(strcmp(info_infix.version, '2.0'));
assert(strcmp(info_infix.isomorphism(1).type, 'equal'));
assert(info_infix.totalIn == 6 && info_infix.totalOut == 4);
assert(isa(X0_infix{1}, 'polytope') && dim(X0_infix{1}) == 6);
assert(isscalar(spec_infix));
assert(dim(spec_infix.set) == 4);
% Y_f[1] < Y_g[1]  →  1 halfspace in the joint output space
assert(size(constraints(spec_infix.set), 1) == 1);

% --- chained equivalence must raise converterIssue ---
try
    vnnlib2cora([CORAROOT '/models/Cora/nn/unitTests/vnnlib/v2/chained_equivalence.vnnlib']);
    assert(false, 'Expected CORAerror for chained equivalence');
catch ME
    assert(contains(ME.identifier, 'CORA'));
    assert(contains(ME.message, 'Chained') || contains(ME.message, 'chained'));
end

res = true;

end

% ------------------------------ END OF CODE ------------------------------
