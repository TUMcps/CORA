function [A,b,Ae,be] = priv_cartProd(A1,b1,Ae1,be1,A2,b2,Ae2,be2)
% priv_cartProd - computes the Cartesian product of two polytopes in
%    halfspace representation
%
% Syntax:
%    [A,b,Ae,be] = priv_cartProd(A1,b1,Ae1,be1,A2,b2,Ae2,be2)
%
% Inputs:
%    A1 - inequality constraint matrix
%    b1 - inequality constraint offset
%    Ae1 - equality constraint matrix
%    be1 - equality constraint offset
%    A2 - inequality constraint matrix
%    b2 - inequality constraint offset
%    Ae2 - equality constraint matrix
%    be2 - equality constraint offset
%
% Outputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%    Ae - equality constraint matrix
%    be - equality constraint offset
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       03-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% block-concatenation of constraint matrices
A = blkdiag(A1,A2);
Ae = blkdiag(Ae1,Ae2);

% vertical concatenation of offset vectors
b = [b1; b2];
be = [be1; be2];

% ------------------------------ END OF CODE ------------------------------
