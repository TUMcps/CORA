function res = testLong_spectraShadow_supportFunc
% testLong_spectraShadow_supportFunc - unit test function of supportFunc
%
% Syntax:
%    res = testLong_spectraShadow_supportFunc
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
% See also: polytope/test_polytope_supportFunc

% Authors:       Adrian Kulmburg
% Written:       14-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty set
SpS = spectraShadow.empty(1);
assert(supportFunc(SpS,1) == -Inf);

% infeasible constraints
A = [-1 -1; -1 1; 1 0];
b = [-2; -2; -1];
P = polytope(A,b);
SpS = spectraShadow(P);

assert(supportFunc(SpS,[1;0]) == -Inf);
assert(supportFunc(SpS,[1;0],'lower') == Inf);


% init unit box as a spectrahedron
n = 3;
A = [eye(n); -eye(n)];
b = [ones(2*n,1)];
P = polytope(A,b);
SpS = spectraShadow(P);

I = eye(n);

for i=1:n
    % direction
    dir = I(:,i);

    % upper bound
    sF_pos = supportFunc(SpS,dir);
    assert(withinTol(sF_pos,1))

    % lower bound
    sF_neg = supportFunc(SpS,dir,'lower');
    assert(withinTol(sF_neg,-1))
end


% init 2D polytope and compute support function along chosen halfspaces
% (note: no redundant halfspaces!)
A = [2 1; 1 3; 2 -2; -2 -2; -1 2];
A = (A' ./ vecnorm(A'))';
b = [2; 1; 2; 2; 2];
P = polytope(A,b);
SpS = spectraShadow(P);

for i=1:length(b)
    sF = supportFunc(SpS,A(i,:)');
    assert(withinTol(sF,b(i)))
end


% unbounded
A = [1 0; -1 0];
b = ones(2,1);
P = polytope(A,b);
SpS = spectraShadow(P);

assert(supportFunc(SpS,[0;1]) == Inf)
assert(supportFunc(SpS,[0;-1]) == Inf)
assert(supportFunc(SpS,[0;1],'lower') == -Inf)
assert(supportFunc(SpS,[0;-1],'lower') == -Inf);


% 2D case with checked support function values
Z = zonotope([1;1],[1.71 -2.14 1.35 0.96; -0.19 -0.84 -1.07 0.12]);
SpS = spectraShadow(Z);

assert(withinTol(supportFunc(SpS,[0;1]),3.22))
assert(withinTol(supportFunc(SpS,[1;0]),7.16))
assert(withinTol(supportFunc(SpS,[0;-1]),1.22))
assert(withinTol(supportFunc(SpS,[-1;0]),5.16))
assert(withinTol(supportFunc(SpS,[1;1]),7.86))
assert(withinTol(supportFunc(SpS,[-1;-1]),3.86))
assert(withinTol(supportFunc(SpS,[-1;1]),6.46))
assert(withinTol(supportFunc(SpS,[1;-1]),6.46))


% polytope without redundant halfspaces
A = [2 1; -2 1; -1 -2];
b = [2; 1; 2];
P = polytope(A,b);
SpS = spectraShadow(P);

% add redundant halfspaces
A_ = [A; [-2 -3; 1 1.5; -2 1.5]];
b_ = [b; 4; 3; 3];
P_ = polytope(A_,b_);
S_ = spectraShadow(P_);

% compare support function evaluations
dirs = [1 1; 4 1; -2 3; -1 -3; 0 -3; -4 2]';
for i=1:size(dirs,2)
    assert(withinTol(supportFunc(SpS,dirs(:,i)),supportFunc(S_,dirs(:,i))))
end


% bounded degenerate polytope (with redundancies)
A = [2 1; -2 1; -1 -2; -2 -3; 2 -1; 1 1.5; -2 1.5];
b = [2; 1; 2; 4; -1; 3; 3];
P = polytope(A,b);
SpS = spectraShadow(P);

% compute support function evaluations
assert(withinTol(supportFunc(SpS,[-2;1],'upper'),supportFunc(SpS,[-2;1],'lower')))
assert(withinTol(supportFunc(SpS,[2;-1],'upper'),supportFunc(SpS,[2;-1],'lower')));


% equality constraints
P = polytope([],[],[1 0],1);
SpS = spectraShadow(P);
assert(withinTol(supportFunc(SpS,[1;0],'upper'),1));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
