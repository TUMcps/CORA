function P_out = plus(P,S)
% plus - overloaded '+' operator for the Minkowski addition of two
%    polytopes or a polytope with a vector
%
% Syntax:
%    P_out = plus(P, S)
%
% Inputs:
%    P - polytope object or numerical vector
%    S - polytope object or numerical vector
%
% Outputs:
%    P_out - polytope after Minkowski addition
%
% Example:
%    A = [2 1; -1 1; -2 -3; 0 -4; 2 -1];
%    b = ones(5,1);
%    P = polytope(A,b);
%    v = [2;1];
%
%    res = P + v;
%
%    figure; hold on
%    plot(P);
%    plot(res,[1,2],'r');
%
% References: MPT-toolbox https://www.mpt3.org/
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Viktor Kotsev
% Written:       20-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% sort arguments (Minkowski addition is commutative)
[P,S] = findClassArg(P,S,'polytope');

% check dimensions
equalDimCheck(P,S);

% dimension
n = dim(P);

%polytope case
if isa(S,"polytope")
    %check for empty polytopes
    if representsa(P, 'emptySet')
        P_out = S; return
    end

    if representsa(S, 'emptySet')
        P_out = P; return
    end

    % all-zero matrices
    PZ = zeros(size(P.A,1),size(S.A,2));
    PZe = zeros(size(P.Ae,1),size(S.Ae,2));

    % inequalities
    A = [S.A -S.A; PZ P.A];
    b = [S.b; P.b];

    % equalities
    Ae = [S.Ae -S.Ae; PZe P.Ae];
    be = [S.be; P.be];

    % project resulting polytope
    temp = polytope(A,b,Ae,be);
    P_out = project(temp,1:n);

elseif isnumeric(S)

    x = S(:);
    A = P.A;
    b = P.b;
    Ae = P.Ae;
    be = P.be;

    if size(A,1) > 0
        b = b + A*x;
    end
    if size(Ae,1) > 0
        be = be + Ae*x;
    end

    P_out = polytope(A,b,Ae,be);
else
    throw(CORAerror('CORA:noops',P,S));

end

% set properties

% If both polytopes are bounded, then sum is also bounded
if (~isempty(P.bounded.val) && P.bounded.val) ...
    && ~isnumeric(S) && (~isempty(S.bounded.val) && S.bounded.val)
    P_out.bounded.val = true;
end

% If one of the polytopes is unbounded, then sum is also unbounded
if (~isempty(P.bounded.val) && ~P.bounded.val) ...
    ||  (~isnumeric(S) && ~isempty(S.bounded.val) && ~S.bounded.val)
    P_out.bounded.val = false;
end

% If one of the polytopes is fully dimensional, then sum is also fully dimensional
if (~isempty(P.fullDim.val) && P.fullDim.val) ...
    ||  (~isnumeric(S) && ~isempty(S.fullDim.val) && S.fullDim.val)
    P_out.fullDim.val = true;
end

end

% ------------------------------ END OF CODE ------------------------------
