function [val,x] = minnorm(Z)
% minnorm - computes the point whose norm is minimal with respect to the
%    center of the zonotope; caution: this function requires the halfspace
%    representation of the zonotope and thus scales exponentially with
%    respect to the dimension
%
% Syntax:
%    [val,x] = minnorm(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    val - value of minimum zonotope norm, i.e., the point on the
%          boundary of Z which has minimum distance to the zonotope center
%    x - point on boundary attaining minimum norm
%
% Example: 
%    Z = zonotope([2;1],[1 -1 0; 1 2 3]);
%    [val,x] = minnorm(Z);
%
% References:
%    [1] S. Boyd and L. Vandenberghe, "Convex optimization", Cambridge
%        University Press, 2004
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: norm

% Authors:       Victor Gassmann
% Written:       18-September-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get halfspace representation
P = polytope(Z - Z.c);
A = P.A; b = P.b;

% compute min norm (obtained by rewriting OP in [1], Sec. 8.4.2, using
% ||a_i||_2 = 1 and argmin -log(det(scalarVar))=argmax scalarVar
[val_squared,ind] = min(b.^2);
x = A(ind,:)'*b(ind) + Z.c;
val = sqrt(val_squared);

% ------------------------------ END OF CODE ------------------------------
