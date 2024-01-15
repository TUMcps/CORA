function P_out = and_(P,S,type,varargin)
% and_ - computes the intersection of a polytope and another set
%    note: the resulting representation is not necessarily minimal!
%
% Syntax:
%    P = and_(P,S,type)
%
% Inputs:
%    P - polytope object
%    S - contSet object or numerical vector
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

% Authors:       Viktor Kotsev, Mark Wetzlinger
% Written:       09-May-2022
% Last update:   14-December-2022 (MW, bug fix, add equality constraints)
%                23-December-2023 (MW, support intersection with numeric)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% re-order such that first argument is polytope
[P,S] = findClassArg(P,S,'polytope');

% quick check: fully empty object -> fullspace
if representsa_(S,'fullspace',0)
    P_out = polytope(P); return
elseif representsa_(P,'fullspace',0)
    P_out = polytope(S); return
end

% check dimension
equalDimCheck(P,S);

% go over cases
switch class(S)
    case 'emptySet'
        % intersection with the empty set yields the empty set
        P_out = polytope.empty(dim(P));

    case 'fullspace'
        % R^n does not impose additional constraints
        P_out = polytope(P);
    
    case 'levelSet'
        % call levelSet function
        P_out = and_(S,P,type);

    case 'double'
        % intersection with numeric: if the point is contained, the
        % intersection is that point, otherwise empty

        % avoid matrices
        if size(S,2) > 1
            throw(CORAerror('CORA:notSupported',...
                'Point clouds not supported for intersection.'));
        end

        if contains_(P,S,'exact',eps)
            P_out = polytope(S);
        else
            P_out = polytope.empty(dim(P));
        end
        
    otherwise
        try
            % convert second object to polytope
            S = polytope(S);
        catch ME
            % no conversion operation implemented
            throw(CORAerror('CORA:noops',P,S));
        end

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

end

% ------------------------------ END OF CODE ------------------------------
