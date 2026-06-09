function res = test_nn_converter_vnnlib2cora_v2()
% test_nn_converter_vnnlib2cora_v2 - native single-network VNN-LIB 2.0
%    features: typed declarations, multi-dimensional indexing, == operator,
%    multiple declare-input blocks per network, safe interval downcast.
%
% Syntax:
%    res = test_nn_converter_vnnlib2cora_v2()
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
% Last update:   08-May-2026 (BK, parser gap coverage: !=, declare-hidden, initialized, bounds checks, equal-to type)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

dir2 = [CORAROOT '/models/Cora/nn/unitTests/vnnlib/v2'];

% --- Single-network 2.0 acasxu_prop_3 ---
[X0, spec, info] = vnnlib2cora([dir2 '/acasxu_prop_3_v2.vnnlib']);
assert(strcmp(info.version, '2.0'));
assert(isscalar(info.networks));
assert(strcmp(info.networks(1).name, 'N'));
assert(isequal(info.networks(1).inputs(1).shape, 5));
assert(isequal(info.networks(1).output.shape,    5));
assert(info.totalIn == 5 && info.totalOut == 5);
assert(isscalar(X0));
assert(isa(X0{1}, 'interval'));
assert(dim(X0{1}) == 5);
assert(isscalar(spec));
assert(strcmp(spec.type, 'unsafeSet'));
assert(all(size(spec.set.A) == [4, 5]));

% --- Multi-dimensional indexing ---
[X0, ~, info] = vnnlib2cora([dir2 '/multidim_basic.vnnlib']);
assert(info.totalIn == 6 && info.totalOut == 2);
assert(isa(X0{1}, 'interval') && dim(X0{1}) == 6);
assert(all(infimum(X0{1}) == 0));
assert(all(supremum(X0{1}) == 1));

% --- Equality operator ---
[X0, ~, info] = vnnlib2cora([dir2 '/equality_singlenet.vnnlib']);
assert(isa(X0{1}, 'interval') && dim(X0{1}) == 3);
inf3 = infimum(X0{1});
sup3 = supremum(X0{1});
assert(inf3(3) == 0.25 && sup3(3) == 0.25, ...
    'equality constraint should pin X[2] to [0.25, 0.25]');

% --- Multiple declare-input blocks (multi-modal pattern) ---
% X1 [2,3]=6 dims + X2 [4]=4 dims -> 10-D joint input, still box -> interval
[X0, ~, info] = vnnlib2cora([dir2 '/multiinput_singlenet.vnnlib']);
assert(strcmp(info.version, '2.0'));
assert(isscalar(info.networks));
assert(numel(info.networks(1).inputs) == 2);
assert(isequal(info.networks(1).inputs(1).shape, [2 3]));
assert(isequal(info.networks(1).inputs(2).shape,  4));
assert(info.totalIn == 10 && info.totalOut == 2);
assert(isscalar(X0));
assert(isa(X0{1}, 'interval'));
assert(dim(X0{1}) == 10);
inf10 = infimum(X0{1});
sup10 = supremum(X0{1});
assert(all(inf10(1:6) == 0)  && all(sup10(1:6) == 1),  'X1 bounds wrong');
assert(all(inf10(7:10) == -1) && all(sup10(7:10) == 1), 'X2 bounds wrong');

% --- Backward-compat 2-output form still works ---
[X0_two, spec_two] = vnnlib2cora([dir2 '/acasxu_prop_3_v2.vnnlib']);
assert(isa(X0_two{1}, 'interval'));
assert(isscalar(spec_two));

% --- Strict < and > operators ---
% strict bounds (0,1)x(0,1) treated as non-strict for polytope representation
[X0, spec] = vnnlib2cora([dir2 '/strict_comparison.vnnlib']);
assert(isa(X0{1}, 'interval') && dim(X0{1}) == 2);
assert(all(infimum(X0{1}) == 0) && all(supremum(X0{1}) == 1));
assert(isscalar(spec));

% --- N-ary arithmetic (+ a b c) and (* scalar var) ---
[X0, spec] = vnnlib2cora([dir2 '/nary_arithmetic.vnnlib']);
assert(isa(X0{1}, 'interval') && dim(X0{1}) == 3);
assert(all(infimum(X0{1}) == 0) && all(supremum(X0{1}) == 1));
% spec should have rows for the n-ary sum and the scaled-variable constraints
assert(isscalar(spec));

% --- Unary negation (- x) ---
% X[0] in [-1,0] expressed via (- X[0]) in [0,1]
[X0, spec] = vnnlib2cora([dir2 '/unary_negate.vnnlib']);
assert(isa(X0{1}, 'interval') && dim(X0{1}) == 2);
infUN = infimum(X0{1}); supUN = supremum(X0{1});
assert(infUN(1) == -1 && supUN(1) == 0, 'unary negate: X[0] should be [-1,0]');
assert(infUN(2) == 0  && supUN(2) == 1, 'unary negate: X[1] should be [0,1]');

% --- Scalar shape [] ---
[X0, spec] = vnnlib2cora([dir2 '/scalar_shape.vnnlib']);
assert(isa(X0{1}, 'interval') && dim(X0{1}) == 1);
assert(infimum(X0{1}) == 0.5 && supremum(X0{1}) == 1.5);

% --- equal-to network equivalence ---
[X0, spec, info] = vnnlib2cora([dir2 '/equal_to.vnnlib']);
assert(strcmp(info.version, '2.0'));
assert(numel(info.networks) == 2);
assert(strcmp(info.networks(2).isomorphicTo, 'f'), 'equal-to should populate isomorphicTo');
assert(strcmp(info.networks(2).equalTo, 'f'), 'equal-to should populate equalTo');
assert(numel(info.isomorphism) == 1);
assert(strcmp(info.isomorphism(1).source, 'g'));
assert(strcmp(info.isomorphism(1).type, 'equal'), 'equal-to should produce type=equal');

% --- != operator: expands to two disjunct input branches ---
[X0, ~] = vnnlib2cora([dir2 '/neq_operator.vnnlib']);
assert(numel(X0) == 2, '!= should produce two disjunct input branches');
assert(isa(X0{1}, 'interval') && dim(X0{1}) == 2);
assert(isa(X0{2}, 'interval') && dim(X0{2}) == 2);
sup1 = supremum(X0{1}); inf2 = infimum(X0{2});
assert(sup1(1) == 0.5, '!= branch 1: upper bound of X[0] should be 0.5');
assert(inf2(1) == 0.5, '!= branch 2: lower bound of X[0] should be 0.5');

% --- declare-hidden: parsed and stored in info.networks.hidden ---
[X0, ~, info] = vnnlib2cora([dir2 '/declare_hidden.vnnlib']);
assert(numel(info.networks(1).hidden) == 1);
assert(strcmp(info.networks(1).hidden(1).name, 'H'));
assert(isequal(info.networks(1).hidden(1).shape, 4));
assert(strcmp(info.networks(1).hidden(1).onnxName, 'hidden_layer_0'));
assert(isa(X0{1}, 'interval') && dim(X0{1}) == 3);

% --- constraint on hidden variable must raise notSupported ---
try
    vnnlib2cora([dir2 '/hidden_constraint.vnnlib']);
    assert(false, 'Expected CORAerror for constraint on hidden variable');
catch ME
    assert(contains(ME.identifier, 'CORA'));
    assert(contains(ME.message, 'hidden') || contains(ME.message, 'H'));
end

% --- index out of bounds must raise converterIssue ---
try
    vnnlib2cora([dir2 '/oob_index.vnnlib']);
    assert(false, 'Expected CORAerror for out-of-bounds index');
catch ME
    assert(contains(ME.identifier, 'CORA'));
    assert(contains(ME.message, 'bounds') || contains(ME.message, 'Index'));
end

% --- duplicate variable declaration must raise converterIssue ---
try
    vnnlib2cora([dir2 '/duplicate_var.vnnlib']);
    assert(false, 'Expected CORAerror for duplicate variable declaration');
catch ME
    assert(contains(ME.identifier, 'CORA'));
    assert(contains(ME.message, 'Multiple declarations') || contains(ME.message, 'X'));
end

% --- initialized flag on declare-input ---
[X0, ~, info] = vnnlib2cora([dir2 '/initialized_input.vnnlib']);
assert(isscalar(info.networks));
assert(numel(info.networks(1).inputs) == 2);
assert(info.networks(1).inputs(1).initialized  == true,  'Xinit should be initialized');
assert(info.networks(1).inputs(2).initialized  == false, 'Xfree should not be initialized');
assert(isa(X0{1}, 'interval') && dim(X0{1}) == 5, 'joint input should be 5-D');

res = true;

end

% ------------------------------ END OF CODE ------------------------------
