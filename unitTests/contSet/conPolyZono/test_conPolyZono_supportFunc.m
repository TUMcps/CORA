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

res = true(0);
tol = 1e-6;

% empty set
cPZ = conPolyZono.empty(2);
res(end+1,1) = supportFunc(cPZ,[1;1],'upper') == -Inf;
res(end+1,1) = supportFunc(cPZ,[1;1],'lower') == Inf;

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
H = 2*blkdiag([temp(2) 0.5*temp(3);0.5*temp(3) temp(4)],0);
f = [temp(1);0;0];
Aeq = [A(2);A(1);A(3)];

options = optimoptions('quadprog','Display','off');
[~,val_] = quadprog(H,f,[],[],Aeq',b,-ones(3,1),ones(3,1),[],options);

% % visualization
% figure; hold on
% plot(cPZ,[1,2],'r','Splits',15);
% ch = conHyperplane(d,val_);
% plot(ch,[1,2],'b');

% compare with exact solution
res(end+1,1) = withinTol(val,val_,tol);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
