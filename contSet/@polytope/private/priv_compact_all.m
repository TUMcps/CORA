function [A,b,Ae,be,empty,minHRep] = priv_compact_all(A,b,Ae,be,n,tol)
% priv_compact_all - removes all redundant constraints in the halfspace
%    representation of an nD polytope
%
% Syntax:
%    [A,b,Ae,be,empty,minHRep] = priv_compact_all(A,b,Ae,be,n,tol)
%
% Inputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%    Ae - equality constraint matrix
%    be - equality constraint offset
%    n - dimension of the polytope
%    tol - tolerance
%
% Outputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%    Ae - equality constraint matrix
%    be - equality constraint offset
%    empty - true/false whether polytope is empty
%    minHRep - minimal representation obtained
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

minHRep = false;

[A,b,Ae,be,empty] = priv_compact_zeros(A,b,Ae,be,tol);
if empty
    return
end

% equality constraints can only be redundant if they are aligned
[A,b,Ae,be] = priv_compact_toEquality(A,b,Ae,be,tol);
[Ae,be,empty] = priv_compact_alignedEq(Ae,be,tol);
if empty
    return
end
[A,b] = priv_compact_alignedIneq(A,b,tol);

% special algorithms for 1D and 2D, general method for nD
if n == 1
    [A,b,Ae,be,empty] = priv_compact_1D(A,b,Ae,be,tol);
elseif n == 2
    [A,b,Ae,be,empty] = priv_compact_2D(A,b,Ae,be,tol);
else
    [A,b,empty] = priv_compact_nD(A,b,Ae,be,n,tol);
end
minHRep = true;

% ------------------------------ END OF CODE ------------------------------
