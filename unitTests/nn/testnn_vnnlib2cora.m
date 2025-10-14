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

% (i) test case: prop_1.vnnlib
% [...] (assert (>= Y_0 3.991125645861615))
prop1Filename = [CORAROOT '/models/Cora/nn/prop_1.vnnlib'];
[X0,specs] = vnnlib2cora(prop1Filename);
% There is only one input set.
assert(isscalar(X0));
% Check input sets constraints.
assert(all(X0{1}.inf == [0.6 -0.5 -0.5 0.45 -0.5]'));
assert(all(X0{1}.sup == [0.679857769 0.5 0.5 0.5 -0.45]'));
% Check output set.
assert(isscalar(specs.set) & isa(specs.set,'polytope') ...
    & representsa(specs.set,'halfspace'));
assert(strcmp(specs.type,'unsafeSet'));
assert(all(withinTol(specs.set.A, [-1 0 0 0 0])));
assert(all(withinTol(specs.set.b, -3.991125645861615)));
% Test samples.
y_u = 2*rand(5,1) - 1; % Random output.
y_u(1) = 3.991125645861615 + abs(y_u(1)); % make it unsafe, i.e., Y_1 >= 3.99...
assert((strcmp(specs.type,'unsafeSet') & contains(specs.set,y_u)) ...
    | (strcmp(specs.type,'safeSet') & ~contains(specs.set,y_u))); % unsafe output.
y_s = 2*rand(5,1) - 1; % Random output.
y_s(1) = 3.991125645861615 - abs(y_s(1)) - 1e-6; % make it safe, i.e., Y_1 < 3.99...
assert((strcmp(specs.type,'unsafeSet') & ~contains(specs.set,y_s)) ...
    | (strcmp(specs.type,'safeSet') & contains(specs.set,y_s))); % safe output.


% (ii) test case: prop_2.vnnlib
% [...] (assert (<= Y_1 Y_0))
%       (assert (<= Y_2 Y_0))
%       (assert (<= Y_3 Y_0))
%       (assert (<= Y_4 Y_0))
prop2Filename = [CORAROOT '/models/Cora/nn/prop_2.vnnlib'];
[X0,specs] = vnnlib2cora(prop2Filename);
% There is only one input set.
assert((isscalar(X0)));
% Check input sets constraints.
assert(all(X0{1}.inf == [0.6 -0.5 -0.5 0.45 -0.5]'));
assert(all(X0{1}.sup == [0.679857769 0.5 0.5 0.5 -0.45]'));
% Check output set.
assert(isscalar(specs.set) & isa(specs.set,'polytope'));
assert(strcmp(specs.type,'unsafeSet'));
assert(all(withinTol(specs.set.A, [-ones(4,1) eye(4)]),'all'));
assert(all(withinTol(specs.set.b, zeros(4,1))));
% Test samples.
y_u = 2*rand(5,1) - 1; % Random output.
y_u(1) = max(y_u(2:5)) + 1e-6; % make it unsafe, i.e., Y_i <= Y_0 for all i=1,...,4
assert((strcmp(specs.type,'unsafeSet') & contains(specs.set,y_u)) ...
    | (strcmp(specs.type,'safeSet') & ~contains(specs.set,y_u))); % unsafe output.
y_s = 2*rand(5,1) - 1; % Random output.
y_s(1) = y_s(randsample(2:5,1)) - 1; % make it safe, i.e., Y_0 < Y_i for any i=1,...,4
assert((strcmp(specs.type,'unsafeSet') & ~contains(specs.set,y_s)) ...
    | (strcmp(specs.type,'safeSet') & contains(specs.set,y_s))); % safe output.

% (iii) test case: prop_3.vnnlib
% [...] (assert (<= Y_0 Y_1))
%       (assert (<= Y_0 Y_2))
%       (assert (<= Y_0 Y_3))
%       (assert (<= Y_0 Y_4))
prop3Filename = [CORAROOT '/models/Cora/nn/unitTests/vnnlib/axas_xu_prop_3.vnnlib'];
[X0,specs] = vnnlib2cora(prop3Filename);
% There is only one input set.
assert(isscalar(X0));
% Check input sets constraints.
assert(all(X0{1}.inf == [-0.303531156 -0.009549297 0.493380324 0.3 0.3]'));
assert(all(X0{1}.sup == [-0.298552812 0.009549297 0.5 0.5 0.5]'));
% Check output set.
assert(isscalar(specs.set) & isa(specs.set,'polytope'));
assert(strcmp(specs.type,'unsafeSet'));
assert(all(withinTol(specs.set.A, [ones(4,1) -eye(4)]),'all'));
assert(all(withinTol(specs.set.b, zeros(4,1))));
% Test samples.
y_u = 2*rand(5,1) - 1; % Random output.
y_u(1) = min(y_u(2:5)) - 1e-6; % make it unsafe, i.e., Y_0 <= Y_i for all i=1,...,4
assert((strcmp(specs.type,'unsafeSet') & contains(specs.set,y_u)) ...
    | (strcmp(specs.type,'safeSet') & ~contains(specs.set,y_u))); % unsafe output.
y_s = 2*rand(5,1) + 1; % Random output.
y_s(1) = y_s(randsample(2:5,1)) + 1e-6; % make it safe, i.e., Y_0 > Y_i for any i=1,...,4
assert((strcmp(specs.type,'unsafeSet') & ~contains(specs.set,y_s)) ...
    | (strcmp(specs.type,'safeSet') & contains(specs.set,y_s))); % safe output.


% (iv) test case: prop_5.vnnlib
% [...] (assert (or
%     (and (<= Y_0 Y_4))
%     (and (<= Y_1 Y_4))
%     (and (<= Y_2 Y_4))
%     (and (<= Y_3 Y_4))
% ))
prop5Filename = [CORAROOT '/models/Cora/nn/prop_5.vnnlib'];
[X0,specs] = vnnlib2cora(prop5Filename);
% There is only one input set.
assert(isscalar(X0));
% Check input sets constraints.
assert(all(X0{1}.inf == [-0.324274257 0.031830989 -0.499999896 -0.5 -0.5]'));
assert(all(X0{1}.sup == [-0.321785085 0.063661977 -0.499204121 -0.227272727 -0.166666667]'));
% Check output set.
assert(isscalar(specs.set) & isa(specs.set,'polytope'));
assert(strcmp(specs.type,'safeSet'));
assert(all(withinTol(specs.set.A, [-eye(4) ones(4,1)]),'all'));
assert(all(withinTol(specs.set.b, zeros(4,1))));
% Test samples.
y_u = 2*rand(5,1) - 1; % Random output.
y_u(5) = y_u(randsample(1:4,1)) + 1e-6; % make it unsafe, i.e., Y_i <= Y_4 for one i=0,...,3
assert(~contains(specs.set,y_u)); % Does not contain the unsafe output.
y_s = 2*rand(5,1) + 1; % Random output.
y_s(5) = min(y_s(1:4)) - 1e-6; % make it safe, i.e., Y_i > Y_4 for all i=0,...,3
assert(contains(specs.set,y_s)); % Does contain the safe output.

% (v) test case: prop_7.vnnlib
% [...] (assert (or
%     (and (<= Y_3 Y_0) (<= Y_3 Y_1) (<= Y_3 Y_2))
%     (and (<= Y_4 Y_0) (<= Y_4 Y_1) (<= Y_4 Y_2))
% ))
prop7Filename = [CORAROOT '/models/Cora/nn/prop_7.vnnlib'];
[X0,specs] = vnnlib2cora(prop7Filename);
% There is only one input set.
assert(isscalar(X0));
% Check input sets constraints.
assert(all(X0{1}.inf == [-0.328422877 -0.499999896 -0.499999896 -0.5 -0.5]'));
assert(all(X0{1}.sup == [0.679857769 0.499999896 0.499999896 0.5 0.5]'));
% Check output set; there are two unsafe sets.
assert(length(specs) == 2);
% Test samples.
for i=1:length(specs)
    % Obtain the i-th specification.
    speci = specs(i);
    assert(isa(speci.set,'polytope'));
    assert(strcmp(speci.type,'unsafeSet'));
    % Create an unsafe output for each of the three cases.
    y_u1 = 2*rand(5,1) - 1; % Random output.
    y_u1(4) = min(y_u1(1:3)) - 1e-6; % make it unsafe (i), i.e., Y_3 <= Y_i for all i=0,...,2
    y_u1(5) = max(y_u1(1:3)) + 1e-6; % ... (not in case ii, i.e., Y_4 > Y_i for any i=0,...,2)
    y_u2 = 2*rand(5,1) - 1; % Random output.
    y_u2(5) = min(y_u2(1:3)) - 1e-6; % make it unsafe (ii), i.e., Y_4 <= Y_i for all i=0,...,2
    y_u2(4) = max(y_u2(1:3)) + 1e-6; % ... (not in case i, i.e., Y_3 > Y_i for any i=0,...,2)
    % Create a safe output.
    y_s1 = 2*rand(5,1) - 1; % Random output.
    y_s1(4) = y_s1(randsample(1:3,1)) + 1e-6; % make it safe, i.e., Y_3 > Y_i and Y_4 > Y_j for any i,j=0,...,2
    y_s1(5) = y_s1(randsample(1:3,1)) + 1e-6; % ...
    % Check the specification: An unsafe set must contain at least one 
    % unsafe output.
    assert(xor(contains(speci.set,y_u1),contains(speci.set,y_u2)));
    % The unsafe set must not contain a safe output.
    assert(~contains(speci.set,y_s1));
end

% (vi) test case: prop_8.vnnlib
% [...] (assert (or
%     (and (<= Y_2 Y_0) (<= Y_2 Y_1))
%     (and (<= Y_3 Y_0) (<= Y_3 Y_1))
%     (and (<= Y_4 Y_0) (<= Y_4 Y_1))
% ))
prop8Filename = [CORAROOT '/models/Cora/nn/prop_8.vnnlib'];
[X0,specs] = vnnlib2cora(prop8Filename);
% There is only one input set.
assert(isscalar(X0));
% Check input sets constraints.
assert(all(X0{1}.inf == [-0.328422877 -0.499999896 -0.015915494 -0.045454545 0.0]'));
assert(all(X0{1}.sup == [0.679857769 -0.374999922 0.015915494 0.5 0.5]'));
% Check output set; there are three unsafe sets.
assert(length(specs) == 3);
% Test samples.
for i=1:length(specs)
    % Obtain the i-th specification.
    speci = specs(i);
    assert(isa(speci.set,'polytope'));
    assert(strcmp(speci.type,'unsafeSet'));
    % Create an unsafe output for each of the three cases.
    y_u1 = 2*rand(5,1) - 1; % Random output.
    y_u1(3) = min(y_u1(1:2)) - 1e-6; % make it unsafe (i), i.e., Y_2 <= Y_i for all i=0,1
    y_u1(4) = max(y_u1(1:2)) + 1e-6; % ... (not in case ii, i.e., Y_3 > Y_i for any i=0,1)
    y_u1(5) = max(y_u1(1:2)) + 1e-6; % ... (not in case iii, i.e., Y_4 > Y_i for any i=0,1)
    y_u2 = 2*rand(5,1) - 1; % Random output.
    y_u2(4) = min(y_u2(1:2)) - 1e-6; % make it unsafe (ii), i.e., Y_3 <= Y_i for all i=0,1
    y_u2(3) = max(y_u2(1:2)) + 1e-6; % ... (not in case i, i.e., Y_2 > Y_i for any i=0,1)
    y_u2(5) = max(y_u2(1:2)) + 1e-6; % ... (not in case iii, i.e., Y_4 > Y_i for any i=0,1)
    y_u3 = 2*rand(5,1) - 1; % Random output.
    y_u3(5) = min(y_u3(1:2)) - 1e-6; % make it unsafe (iii), i.e., Y_4 <= Y_i for all i=0,1
    y_u3(3) = max(y_u3(1:2)) + 1e-6; % ... (not in case i, i.e., Y_2 > Y_i for any i=0,1)
    y_u3(4) = max(y_u3(1:2)) + 1e-6; % ... (not in case ii, i.e., Y_3 > Y_i for any i=0,1)
    % Create a safe output.
    y_s1 = 2*rand(5,1) - 1; % Random output.
    y_s1(3) = y_s(randsample(1:2,1)) + 1e-6; % make it safe, i.e., Y_2 > Y_i, Y_3 > Y_j, and Y_4 > Y_k for any i,j,k=0,1
    y_s1(4) = y_s(randsample(1:2,1)) + 1e-6; % ...
    y_s1(5) = y_s(randsample(1:2,1)) + 1e-6; % ...
    % Check if the specification: An unsafe set must contain at least one 
    % unsafe output.
    assert(xor(contains(speci.set,y_u1), ...
        xor(contains(speci.set,y_u2),contains(speci.set,y_u3))));
    % An unsafe set must not contain a safe output.
    assert(~contains(speci.set,y_s1));
end

% (vii) test case: prop_9.vnnlib
% [...] (assert (or
%     (and (<= Y_0 Y_3))
%     (and (<= Y_1 Y_3))
%     (and (<= Y_2 Y_3))
%     (and (<= Y_4 Y_3))
% ))
prop9Filename = [CORAROOT '/models/Cora/nn/prop_9.vnnlib'];
[X0,specs] = vnnlib2cora(prop9Filename);
% There is only one input set.
assert(isscalar(X0));
% Check input sets constraints.
assert(all(X0{1}.inf == [-0.295233916 -0.063661977 -0.499999896 -0.5 -0.5]'));
assert(all(X0{1}.sup == [-0.212261512 -0.022281692 -0.498408347 -0.454545455 -0.375]'));
% Check output set.
assert(isscalar(specs.set) & (isa(specs.set,'polytope')));
assert(strcmp(specs.type,'safeSet'));
A = [-eye(4), ones(4,1)];
assert(all(withinTol(specs.set.A, A(:,[1:3 5 4])),'all'));
assert(all(withinTol(specs.set.b, zeros(4,1))));
% Test samples.
y_u = 2*rand(5,1) - 1; % Random output.
y_u(4) = y_u(randsample([1:3 5],1)) + 1; % make it unsafe
assert((strcmp(specs.type,'unsafeSet') & contains(specs.set,y_u)) ...
    | (strcmp(specs.type,'safeSet') & ~contains(specs.set,y_u))); % unsafe output.
y_s = 2*rand(5,1) + 1; % Random output.
y_s(4) = min(y_s([1:3 5])) - 1; % make it safe
assert((strcmp(specs.type,'unsafeSet') & ~contains(specs.set,y_s)) ...
    | (strcmp(specs.type,'safeSet') & contains(specs.set,y_s))); % safe output.

end

% ------------------------------ END OF CODE ------------------------------
