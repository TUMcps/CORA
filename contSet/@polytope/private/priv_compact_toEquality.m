function [A,b,Ae,be] = priv_compact_toEquality(A,b,Ae,be,tol)
% priv_compact_toEquality - rewrite all pairwise inequality constraints
%    Ax <= b, -Ax <= -b as equality constraints Ax == b
%    note: expects normalized constraints with respect to 'A'
%
% Syntax:
%    [A,b,Ae,be] = priv_compact_toEquality(A,b,Ae,be,tol)
%
% Inputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%    Ae - equality constraint matrix
%    be - equality constraint offset
%    tol - tolerance
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

% concatenate normal vectors and offsets
Ab = [A, b];

nrIneq = size(Ab,1);
idxMoveToEq = false(1,nrIneq);
for i=1:nrIneq
    % skip constraint if already matched
    if ~idxMoveToEq(i)

        % check if there exists a constraint with an inverted normal vector
        % by computing the sum of the entries (must be zero up to tolerance)
        sum_of_vectors = Ab' + Ab(i,:)';
        sum_of_vectors(withinTol(sum_of_vectors,0,tol)) = 0;
        idxInverted = ~any(sum_of_vectors,1);
        
        % search for a match
        if any(idxInverted)
            idxMoveToEq = idxMoveToEq | idxInverted;
            idxMoveToEq(i) = true;
            Ae = [Ae; A(i,:)];
            be = [be; b(i)];
        end
    end
end

% remove all pairwise inequality constraints
A = A(~idxMoveToEq,:);
b = b(~idxMoveToEq);

% ------------------------------ END OF CODE ------------------------------
