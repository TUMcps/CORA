function res = testnn_vnnlib2cora()
% testnn_vnnlib2cora - tests reading vnnlib files
%
% Syntax:
%    res = testnn_vnnlib2cora()
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
% See also: -

% Authors:       Lukas Koller
% Written:       03-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% We use the specs from the acasxu benchmark: prop_1, prop_2, prop_3, and
% prop_5.

res = true;

% First test case: prop_1.vnnlib
prop1Filename = [CORAROOT '/models/Cora/nn/prop_1.vnnlib'];
[X0,specs] = vnnlib2cora(prop1Filename);
% There is only one input set.
assert((isscalar(X0)));
% Check input sets constraints.
assert(all(X0{1}.inf == [0.6 -0.5 -0.5 0.45 -0.5]'));
assert(all(X0{1}.sup == [0.679857769 0.5 0.5 0.5 -0.45]'));
% Check output set.
assert((isscalar(specs.set)) & isa(specs.set,'polytope') & representsa(specs.set,'halfspace'));
assert((strcmp(specs.type,'unsafeSet')));
assert(all(withinTol(specs.set.A, [-1 0 0 0 0])));
assert(all(withinTol(specs.set.b, -3.991125645861615)));


% Second test case: prop_2.vnnlib
prop2Filename = [CORAROOT '/models/Cora/nn/prop_2.vnnlib'];
[X0,specs] = vnnlib2cora(prop2Filename);
% There is only one input set.
assert((isscalar(X0)));
% Check input sets constraints.
assert(all(X0{1}.inf == [0.6 -0.5 -0.5 0.45 -0.5]'));
assert(all(X0{1}.sup == [0.679857769 0.5 0.5 0.5 -0.45]'));
% Check output set.
assert((isscalar(specs.set)) & (isa(specs.set,'polytope')));
assert((strcmp(specs.type,'unsafeSet')));
assert(all(withinTol(specs.set.A, [-ones(4,1) eye(4)]),'all'));
assert(all(withinTol(specs.set.b, zeros(4,1))));


% Third test case: prop_3.vnnlib
prop3Filename = [CORAROOT '/models/Cora/nn/unitTests/vnnlib/axas_xu_prop_3.vnnlib'];
[X0,specs] = vnnlib2cora(prop3Filename);
% There is only one input set.
assert((isscalar(X0)));
% Check input sets constraints.
assert(all(X0{1}.inf == [-0.303531156 -0.009549297 0.493380324 0.3 0.3]'));
assert(all(X0{1}.sup == [-0.298552812 0.009549297 0.5 0.5 0.5]'));
% Check output set.
assert((isscalar(specs.set)) & (isa(specs.set,'polytope')));
assert((strcmp(specs.type,'unsafeSet')));
assert(all(withinTol(specs.set.A, [ones(4,1) -eye(4)]),'all'));
assert(all(withinTol(specs.set.b, zeros(4,1))));


% Fourth test case: prop_5.vnnlib
prop5Filename = [CORAROOT '/models/Cora/nn/prop_5.vnnlib'];
[X0,specs] = vnnlib2cora(prop5Filename);
% There is only one input set.
assert((isscalar(X0)));
% Check input sets constraints.
assert(all(X0{1}.inf == [-0.324274257 0.031830989 -0.499999896 -0.5 -0.5]'));
assert(all(X0{1}.sup == [-0.321785085 0.063661977 -0.499204121 -0.227272727 -0.166666667]'));
% Check output set.
assert((isscalar(specs.set)) & (isa(specs.set,'polytope')));
assert((strcmp(specs.type,'safeSet')));
assert(all(withinTol(specs.set.A, [-eye(4) ones(4,1)]),'all'));
assert(all(withinTol(specs.set.b, zeros(4,1))));

end

% ------------------------------ END OF CODE ------------------------------
