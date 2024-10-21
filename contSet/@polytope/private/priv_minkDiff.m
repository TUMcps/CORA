function [A,b,Ae,be,empty] = priv_minkDiff(A,b,Ae,be,S)
% priv_minkDiff - computes the Minkowski difference of a polytope and
%    another set or point using support functions
%
% Syntax:
%    [A,b,Ae,be,empty] = priv_minkDiff(A,b,Ae,be,S)
%
% Inputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%    Ae - equality constraint matrix
%    be - equality constraint offset
%    S - contSet object
%
% Outputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%    Ae - equality constraint matrix
%    be - equality constraint offset
%    empty - true/false whether result is the empty set
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

% assume non-empty result
empty = false;

% shift entry in offset vector by support function value of subtrahend
for i = 1:size(A,1)
    l = supportFunc_(S,A(i,:)','upper');
    if isinf(l)
        % subtrahend is unbounded in a direction where the minuend is
        % bounded -> result is empty
        A = []; b = []; Ae = []; be = [];
        empty = true;
        return
    end
    b(i) = b(i) - l;
end

% for equality constraints, we can easily detect if the resulting set
% becomes empty, as the lower and upper bounds for the support function
% in those directions must be equal
for i = 1:size(Ae,1)
    I = supportFunc_(S,Ae(i,:)','range');
    if ~withinTol(I.inf,I.sup)
        A = []; b = []; Ae = []; be = [];
        empty = true;
        return
    end
    be(i) = be(i) - I.inf;
end

% ------------------------------ END OF CODE ------------------------------
