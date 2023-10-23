function res = test_conZonotope_cartProd
% test_conZonotope_cartProd - unit test function for Cartesian product
%
% Syntax:
%    res = test_conZonotope_cartProd
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
% Written:       03-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% constrained zonotope - constrained zonotope case:

% constrained zonotope 1
c1 = [0; 0];
G1 = [3 0 1; 0 2 1];
A1 = [1 0 1]; b1 = 1;
cZ1 = conZonotope(c1,G1,A1,b1);

% constrained zonotope 2
c2 = [1; 0];
G2 = [1.5 -1.5 0.5; 1 0.5 -1];
A2 = [1 1 1]; b2 = 1;
cZ2 = conZonotope(c2,G2,A2,b2);

% compute Cartesian product
cZ_ = cartProd(cZ1,cZ2);

% true result
cZ_true = conZonotope([c1;c2],blkdiag(G1,G2),blkdiag(A1,A2),[b1;b2]);

% compare results
if ~isequal(cZ_,cZ_true)
    res = false;
end


% constrained zonotope - zonotope case:

% zonotope 2
c2 = [1; 0];
G2 = [1.5 -1.5 0.5; 1 0.5 -1];
Z2 = zonotope(c2,G2);

% compute Cartesian product
cZ_ = cartProd(cZ1,Z2);

% true result
cZ_true = conZonotope([c1;c2],blkdiag(G1,G2),[A1,zeros(1,3)],b1);

% compare results
if ~isequal(cZ_,cZ_true)
    res = false;
end


% constrained zonotope - interval case:

% interval 2
lb2 = [2; 3];
ub2 = [4; 3];
I2 = interval(lb2,ub2);

% compute Cartesian product
cZ_ = cartProd(cZ1,I2);

% true result
cZ_true = conZonotope([c1;0.5*(lb2+ub2)],blkdiag(G1,0.5*diag(ub2-lb2)),...
    [A1,zeros(1,2)],b1);

% compare results
if ~isequal(cZ_,cZ_true)
    res = false;
end


% constrained zonotope - numeric case (two orderings):

% numeric 2
num = 1;

% compute Cartesian product
cZ_ = cartProd(cZ1,num);
cZ__ = cartProd(num,cZ1);

% true results
cZ_true = conZonotope([c1;num],[G1;zeros(1,3)],A1,b1);
cZ__true = conZonotope([num;c2],[zeros(1,3);G1],A1,b1);

% compare results
if ~isequal(cZ_,cZ_true) || ~isequal(cZ__,cZ__true)
    res = false;
end


% ------------------------------ END OF CODE ------------------------------
