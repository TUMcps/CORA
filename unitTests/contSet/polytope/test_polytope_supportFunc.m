function res = test_polytope_supportFunc
% test_polytope_supportFunc - unit test function of supportFunc
%
% Syntax:
%    res = test_polytope_supportFunc
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
% See also: none

% Authors:       Mark Wetzlinger
% Written:       29-November-2022
% Last update:   15-November-2023 (MW, vertex representation case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1D, bounded
A = [1;-1]; b = [2;0.5];
P = polytope(A,b);
res(end+1,1) = supportFunc(P,1) == 2;
res(end+1,1) = supportFunc(P,1,'lower') == -0.5;

% 1D, unbounded
A = 1; b = 1;
P = polytope(A,b);
res(end+1,1) = supportFunc(P,1) == 1;
res(end+1,1) = supportFunc(P,-1) == Inf;

% 1D, fully empty
A = zeros(0,1); b = zeros(0,0);
P = polytope(A,b);
res(end+1,1) = supportFunc(P,1) == Inf;
res(end+1,1) = supportFunc(P,-1) == Inf;


% 2D, infeasible constraints -> empty
A = [-1 -1; -1 1; 1 0]; b = [-2; -2; -1];
P = polytope(A,b);
res(end+1,1) = supportFunc(P,[1;0]) == -Inf;
res(end+1,1) = supportFunc(P,[1;0],'lower') == Inf;

% 2D, vertex representation
V = [1 0; 0.8 1; 0 1.4; -0.5 0.5; 0.2 -1]';
P = polytope(V);
[val,x] = supportFunc(P,[1;0]);
res(end+1,1) = withinTol(val,1) && all(withinTol(x,[1;0]));
[val,x] = supportFunc(P,[1;0],'lower');
res(end+1,1) = withinTol(val,-0.5) && all(withinTol(x,[-0.5;0.5]));
[val,x] = supportFunc(P,[1;0],'range');
res(end+1,1) = isequal(val,interval(-0.5,1)) ...
    && all(all(withinTol(x,[-0.5 0.5; 1 0]')));

% 2D, no redundant halfspaces
A = [2 1; 1 3; 2 -2; -2 -2; -1 2];
A = (A' ./ vecnorm(A'))';
b = [2; 1; 2; 2; 2];
P = polytope(A,b);
% compute support function along chosen halfspaces
for i=1:length(b)
    sF = supportFunc(P,A(i,:)');
    if ~withinTol(sF,b(i))
        throw(CORAerror('CORA:testFailed'));
    end
end

% 2D, unbounded
A = [1 0; -1 0]; b = ones(2,1);
P = polytope(A,b);
res(end+1,1) = supportFunc(P,[0;1]) == Inf ...
        && supportFunc(P,[0;-1]) == Inf ...
        && supportFunc(P,[0;1],'lower') == -Inf ...
        && supportFunc(P,[0;-1],'lower') == -Inf;

% 2D, comparison to support function values from other representation
Z = zonotope([1;1],[1.71 -2.14 1.35 0.96; -0.19 -0.84 -1.07 0.12]);
P = polytope(Z);
res(end+1,1) = withinTol(supportFunc(P,[0;1]),3.22) ...
        && withinTol(supportFunc(P,[1;0]),7.16) ...
        && withinTol(supportFunc(P,[0;-1]),1.22) ...
        && withinTol(supportFunc(P,[-1;0]),5.16) ...
        && withinTol(supportFunc(P,[1;1]),7.86) ...
        && withinTol(supportFunc(P,[-1;-1]),3.86) ...
        && withinTol(supportFunc(P,[-1;1]),6.46) ...
        && withinTol(supportFunc(P,[1;-1]),6.46);

% 2D, without/with redundant halfspaces
A = [2 1; -2 1; -1 -2];
b = [2; 1; 2];
P = polytope(A,b);
A_ = [A; [-2 -3; 1 1.5; -2 1.5]];
b_ = [b; 4; 3; 3];
P_ = polytope(A_,b_);
% compare support function evaluations
dirs = [1 1; 4 1; -2 3; -1 -3; 0 -3; -4 2]';
for i=1:size(dirs,2)
    if ~withinTol(supportFunc(P,dirs(:,i)),supportFunc(P_,dirs(:,i)))
        throw(CORAerror('CORA:testFailed'));
    end
end

% 2D, bounded, degenerate, with redundancies
A = [2 1; -2 1; -1 -2; -2 -3; 2 -1; 1 1.5; -2 1.5];
b = [2; 1; 2; 4; -1; 3; 3];
P = polytope(A,b);
% compute support function evaluations
res(end+1,1) = withinTol(supportFunc(P,[-2;1],'upper'),supportFunc(P,[-2;1],'lower')) ...
        && withinTol(supportFunc(P,[2;-1],'upper'),supportFunc(P,[2;-1],'lower'));

% 2D, only equality constraints
Ae = [1 0]; be = 1;
P = polytope([],[],Ae,be);
res(end+1,1) = withinTol(supportFunc(P,[1;0],'upper'),1);

% 2D, trivially fulfilled constraints
A = [0 0]; b = 1; Ae = [0 0]; be = 0;
P = polytope(A,b,Ae,be);
res(end+1,1) = supportFunc(P,[1;1]) == Inf;
res(end+1,1) = supportFunc(P,[1;1],'lower') == -Inf;


% 3D, box
n = 3; A = [eye(n); -eye(n)]; b = [ones(2*n,1)];
P = polytope(A,b);
for i=1:n
    % direction
    ei = zeros(n,1); ei(i) = 1;
    % upper bound and lower bound
    if ~withinTol(supportFunc(P,ei),1) || ~withinTol(supportFunc(P,ei,'lower'),-1)
        throw(CORAerror('CORA:testFailed'));
    end
end


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
