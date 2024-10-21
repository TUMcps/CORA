function res = test_conPolyZono_supportFunc
% test_conPolyZono_supportFunc - unit test function for support function 
%                                enclosure of constrained polynomial 
%                                zonotopes
%
% Syntax:
%    test_conPolyZono_supportFunc
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
% See also: conPolyZono/supportFunc

% Authors:       Niklas Kochdumper
% Written:       26-January-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-6;

% empty set
cPZ = conPolyZono.empty(2);
assert(supportFunc(cPZ,[1;1],'upper') == -Inf);
assert(supportFunc(cPZ,[1;1],'lower') == Inf);

% Analytical Tests --------------------------------------------------------

% define constrained polynomial zonotope
c = [0;0];
G = [1 0 1 -1; 0 2 1 2];
E = [1 2 1 0; 0 0 1 2; 0 0 0 0];
A = [1 1 0.5];
b = 0.5;
EC = [0 1 0;1 0 0; 0 0 1];

cPZ = conPolyZono(c,G,E,A,b,EC);

% define direction
d = [1;1];

% remove constraints
val = supportFunc(cPZ,d,'lower','quadProg');

% compute exact solution
temp = d'*G;
problem.H = 2*blkdiag([temp(2) 0.5*temp(3);0.5*temp(3) temp(4)],0);
problem.f = [temp(1);0;0];
problem.Aineq = [];
problem.bineq = [];
problem.Aeq = [A(2);A(1);A(3)]';
problem.beq = b;
problem.lb = -ones(3,1);
problem.ub = ones(3,1);

[~,val_] = CORAquadprog(problem);

% % visualization
% figure; hold on
% plot(cPZ,[1,2],'r','Splits',15);
% P = polytope([],[],d,val_);
% plot(P,[1,2],'b');

% compare with exact solution
assert(withinTol(val,val_,tol));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
