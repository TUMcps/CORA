function P_out = minkDiff(P,S,varargin)
% minkDiff - compute the Minkowski difference of two polytopes:
%         P - S = P_out <-> P_out + S \subseteq P
%
% Syntax:
%    P_out = minkDiff(P,S)
%    P_out = minkDiff(P,S,type)
%
% Inputs:
%    P - polytope object
%    S - contSet object, or numerical vector
%    type - type of computation
%           'exact': using support function evaluation along directions of
%                    the halfspaces of the minuend
%           'exact:vertices': using vertices and intersection
%
% Outputs:
%    P_out - polytope object after Minkowski difference
%
% Example: 
%    P1 = polytope([1 0 -1 0 1;0 1 0 -1 1]',[4;4;4;4;4]);
%    P2 = polytope([-1 0 1;-2 3 0]);
%
%    P = minkDiff(P1,P2);
%
%    figure; hold on;
%    plot(P1);
%    plot(P2,[1,2],'r');
%    plot(P,[1,2],'g');
%
% References:
%    [1] M. Althoff, "On Computing the Minkowski Difference of Zonotopes"
%    [2] M. Wetzlinger, V. Kotsev, A. Kulmburg, M. Althoff. "Implementation
%        of Polyhedral Operations in CORA 2024", ARCH'24.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/minkDiff

% Authors:       Niklas Kochdumper
% Written:       04-February-2021
% Last update:   01-December-2022 (MW, support function method as default)
%                23-November-2023 (MW, bug fix for equality constraints)
%                14-July-2024 (MW, bug fix for equality constraints)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set default values
type = setDefaultValues({'exact'},varargin);

% check input arguments
inputArgsCheck({{P,'att','polytope'}, ...
                {S,'att',{'contSet','numeric'}},...
                {type,'str',{'exact','exact:vertices'}}});

% different algorithms for different set representations
if isnumeric(S)
    P_out = P + (-S);
    return    
end

% check dimensions
equalDimCheck(P,S);

% read out dimension
n = dim(P);

% fullspace minuend
if representsa_(P,'fullspace',0)
    P_out = polytope.Inf(n);
    return
end

% 1D case for polytopes
if n == 1 && isa(S,'polytope')
    P_out = aux_minkDiff_1D(P,S);
    return
end

% exact computation
if strcmp(type,'exact')
    % only requires support function evaluation of subtrahend, see [1, (27)]

    % ensure that halfspace representation is there
    constraints(P);

    % read out constraints
    A = P.A_.val; b = P.b_.val;
    Ae = P.Ae_.val; be = P.be_.val;
    
    % shift entry in offset vector by support function value of subtrahend
    for i = 1:size(A,1)
        l = supportFunc_(S,A(i,:)','upper');
        if isinf(l)
            % subtrahend is unbounded in a direction where the minuend is
            % bounded -> result is empty
            P_out = polytope.empty(n); return
        end
        b(i) = b(i) - l;
    end
    % for equality constraints, we can easily detect if the resulting set
    % becomes empty, as the lower and upper bounds for the support function
    % in those directions must be equal
    for i = 1:size(Ae,1)
        I = supportFunc_(S,Ae(i,:)','range');
        if ~withinTol(I.inf,I.sup)
            P_out = polytope.empty(n); return
        end
        be(i) = be(i) - I.inf;
    end
    
    % init resulting polytope
    P_out = polytope(A,b,Ae,be);

elseif strcmp(type,'exact:vertices')
    if isa(S,'zonotope') || isa(S,'interval')
    
        if isa(S,'interval')
            % convert to zonotope
            Z = zonotope(S);
        end
        
        % compute Minkowski diff. according to Theorem 1 in [1]
        c = center(Z);
        G = generators(Z);
        
        % we will require the H-representation for the intersections below
        P_out = P1 + (-c);
        constraints(P_out);
        
        for i = 1:size(G,2)
            P_out = and_(P_out + G(:,i),P_out + (-G(:,i)),'exact');
        end
    
    elseif isa(S,'conZonotope') || isa(S,'polytope') || ...
            isa(S,'zonoBundle')
    
        % compute Minkowski difference according to Lemma 1 in [1]
        V = vertices(S);
        P_out = P + (-V(:,1));
        constraints(P_out);
        
        for i = 2:size(V,2)
            P_out = and_(P_out, P + (-V(:,i)),'exact'); 
        end
     
    else

        throw(CORAerror('CORA:noops',P,S));
    end
end

% set properties
P_out = aux_setproperties(P_out,P,S);

end


% Auxiliary functions -----------------------------------------------------

function P_out = aux_minkDiff_1D(P,S)

    % dimension
    n = dim(P);

    % remove redundancies -> we get one of three cases:
    % - two inequality constraints
    % - empty set
    % - one equality constraint
    P_ = compact_(P,'all',1e-9);
    % if the minuend is fullspace, minkDiff(P1,S) = fullspace
    if representsa_(P_,'fullspace',0)
        P_out = polytope.Inf(n); return
    end

    S_ = compact(S);
    % if the subtrahend is empty, minkDiff(P1,S) = R^n
    if representsa_(S_,'emptySet',1e-10)
        P_out = fullspace(n); return
    end

    if ~isempty(P_.Ae_.val) || ~isempty(S_.Ae_.val)
        % currently not supported -> fix
        throw(CORAerror('CORA:notSupported',...
            'minkDiff for equality constraints currently not supported.'));
    end

    % both A matrices are normalized in minHRep (only 1D)
    P_out = polytope(P_.A_.val, P_.b_.val - S_.b_.val);

end

function P_out = aux_setproperties(P_out,P,S)

% currently only polytopes supported
if isa(S,'polytope')

    % If both polytopes are bounded, then difference is also bounded
    if (~isempty(P.bounded.val) && P.bounded.val) ...
        && ~isnumeric(S) && (~isempty(S.bounded.val) && S.bounded.val)
        P_out.bounded.val = true;
    end
    
    % If one of the polytopes is unbounded, then difference is also unbounded
    if (~isempty(P.bounded.val) && ~P.bounded.val) ...
        || (~isnumeric(S) && ~isempty(S.bounded.val) && ~S.bounded.val)
        P_out.bounded.val = false;
    end
    
    % If one of the polytopes is fully dimensional, then difference is also fully dimensional
    if (~isempty(P.fullDim.val) && P.fullDim.val) ...
        ||  (~isnumeric(S) && ~isempty(S.fullDim.val) && S.fullDim.val)
        P_out.fullDim.val = true;
    end

end

end

% ------------------------------ END OF CODE ------------------------------
