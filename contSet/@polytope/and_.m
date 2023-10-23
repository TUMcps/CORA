function P_out = and_(P,S,type,varargin)
% and_ - computes the intersection of a polytope and another set
%    note: the resulting representation is not necessarily minimal!
%
% Syntax:
%    P = and_(P,S,type)
%
% Inputs:
%    P - polytope object
%    S - contSet object
%    type - 'exact' or 'approx'
%
% Outputs:
%    P - polytope object
%
% Example: 
%    P = polytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    hyp = conHyperplane([1 1],2,[-1 0],-1);
%
%    res = P & hyp;
%
%    figure; hold on
%    xlim([-2,4]); ylim([-4,4]);
%    plot(hyp,[1,2],'r','LineWidth',3);
%    plot(P,[1,2],'b');
%    plot(res,[1,2],'g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/and, conZonotope/and_

% Authors:       Viktor Kotsev
% Written:       09-May-2022
% Last update:   14-December-2022 (MW, bug fix, add equality constraints)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% re-order such that first argument is polytope
[P,S] = findClassArg(P,S,'polytope');

% quick check: empty set case
if isemptyobject(S) || isemptyobject(P)
    P_out = polytope(zeros(0,dim(P)),[]);
    % set properties set in constructor
    return;
end

% check dimension
equalDimCheck(P,S);

% check fullspace
if isa(S,'fullspace')
    P_out = polytope(P);
    return;
end

% call levelSet method
if isa(S,'levelSet')
    P_out = and_(S,P,type);
    return;
end

% convert second object to polytope
S = polytope(S);

% compute intersection
P_out = polytope([P.A; S.A], [P.b; S.b], [P.Ae; S.Ae], [P.be; S.be]);

% set properties
% intersection with empty set is empty
if (~isempty(P.emptySet.val) && P.emptySet.val) ...
        || (~isempty(S.emptySet.val) && S.emptySet.val)
    P_out.emptySet.val = true;
end

% intersection with bounded set yields a bounded set
if (~isempty(P.bounded.val) && P.bounded.val) ...
        || (~isempty(S.bounded.val) && S.bounded.val)
    P_out.bounded.val = true;
end

% intersection with a degenerate set yields a degenerate set
if (~isempty(P.fullDim.val) && ~P.fullDim.val) ...
        || (~isempty(S.fullDim.val) && ~S.fullDim.val)
    P_out.fullDim.val = false;
end

% ------------------------------ END OF CODE ------------------------------
