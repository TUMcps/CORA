function P_out = cartProd_(P,S,type,varargin)
% cartProd_ - returns the Cartesian product of a polytope and another set
%    or point
% 
% Syntax:
%    P_out = cartProd_(P,S,type,varargin)
%
% Inputs:
%    P - polytope object, numerical vector
%    S - contSet object, numerical vector
%    type - 'exact', 'inner', 'outer'
%
% Outputs:
%    P_out - resulting polytope object
%
% Example: 
%    P1 = polytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    P2 = polytope([1;-1],[3;2]);
%    P = cartProd(P1,P2);
%
% Reference:
%    [1] M. Wetzlinger, V. Kotsev, A. Kulmburg, M. Althoff. "Implementation
%        of Polyhedral Operations in CORA 2024", ARCH'24.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/cartProd, zonotope/cartProd_

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       26-November-2019
% Last update:   05-May-2020 (MW, standardized error message)
%                04-April-2022 (VK, adapted to polytope class)
% Last revision: 14-July-2024 (MW, refactor, support V-Rep)

% ------------------------------ BEGIN CODE -------------------------------

% ensure that not a Cartesian product with a point cloud (should be fast)
if (isnumeric(P) && size(P,2) > 1) || (isnumeric(S) && size(S,2) > 1)
    throw(CORAerror('CORA:notSupported',...
        'Cartesian product only supported for single points.'));
end

% note: no re-ordering!

% first or second set is polytope
if isa(P,'polytope')

    % different cases for different set representations
    if isa(S,'polytope')
        % polytope x polytope
        if P.isHRep.val 
            if S.isHRep.val
                P_out = aux_cartProd_Hpoly_Hpoly(P,S);
            elseif S.isVRep.val
                P_out = aux_cartProd_Hpoly_Vpoly(P,S,true);
            end
        elseif P.isVRep.val
            if S.isHRep.val
                P_out = aux_cartProd_Hpoly_Vpoly(S,P,false);
            elseif S.isVRep.val
                P_out = aux_cartProd_Vpoly_Vpoly(P,S);
            end
        end
        % infer set properties
        P_out = aux_cartProd_polytope_properties(P_out,P,S);
        
    elseif isnumeric(S)
        % polytope x point
        P_out = aux_cartProd_polytope_point(P,S,true);

    elseif isa(S,'zonoBundle') || isa(S,'zonotope') || ...
           isa(S,'interval') || isa(S,'conZonotope')
        P_out = cartProd_(P,polytope(S),type);
        
    elseif isa(S,'polyZonotope')
        P_out = cartProd_(polyZonotope(P),S,type);
        
    else
        % throw error for given arguments
        throw(CORAerror('CORAerror:noops',P,S));
    end

else

    % different cases for different set representations
    if isnumeric(P)
        % point x polytope (... second set has to be a polytope, otherwise
        % this function would not have been called)
        P_out = aux_cartProd_polytope_point(S,P,false);

    else
        % throw error for given arguments
        throw(CORAerror('CORAerror:noops',P,S));
    end  
end

end


% Auxiliary functions -----------------------------------------------------

function P_out = aux_cartProd_Hpoly_Hpoly(P,S)
% Cartesian product of two polytopes in halfspace representation according 
% to [1, (29)]

% block-concatenation of constraint matrices
A = blkdiag(P.A_.val,S.A_.val);
Ae = blkdiag(P.Ae_.val,S.Ae_.val);
% vertical concatenation of offset vectors
b = [P.b_.val; S.b_.val];
be = [P.be_.val; S.be_.val];

% init resulting polytope
P_out = polytope(A,b,Ae,be);

end

function P_out = aux_cartProd_Hpoly_Vpoly(P_H,P_V,order)
% Cartesian product of a polytope P_H in halfspace representation and a
% polytope P_V in vertex representation
% helper variable: order = true  (H-poly x V-poly)
%                  order = false (V-poly x H_poly)

% we convert the vertex representation into the halfspace representation
constraints(P_V);
if order
    P_out = aux_cartProd_Hpoly_Hpoly(P_H,P_V);
else
    P_out = aux_cartProd_Hpoly_Hpoly(P_V,P_H);
end

end

function P_out = aux_cartProd_Vpoly_Vpoly(P1,P2)
% Cartesian product of two polytopes in vertex representation according to
% [1, (28)]

% read out number of vertices
m1 = size(P1.V_.val,2);
m2 = size(P2.V_.val,2);

% concatenate vertex matrices
V_cartProd = [repelem(P1.V_.val,1,m2); repmat(P2.V_.val,1,m1)];

% instantiate resulting Cartesian product (also in V-representation)
P_out = polytope(V_cartProd);

end

function P_out = aux_cartProd_polytope_properties(P_out,P1,P2)
% infer properties of Cartesian product according to [1, Table 1]; as a
% rule of thumb, the 'more restrictive' property wins, e.g., empty x
% non-empty/unknown results in empty
% note: order of polytopes (or representation) does not matter here

% emptiness
if (~isempty(P1.emptySet.val) && P1.emptySet.val) ...
        || (~isempty(P2.emptySet.val) && P2.emptySet.val)
    P_out.emptySet.val = true;
elseif (~isempty(P1.emptySet.val) && ~P1.emptySet.val) ...
        && (~isempty(P2.emptySet.val) && ~P2.emptySet.val)
    P_out.emptySet.val = false;
end

% boundedness
if (~isempty(P1.bounded.val) && P1.bounded.val) ...
        || (~isempty(P2.bounded.val) && P2.bounded.val)
    P_out.bounded.val = true;
elseif (~isempty(P1.bounded.val) && ~P1.bounded.val) ...
        && (~isempty(P2.bounded.val) && ~P2.bounded.val)
    P_out.bounded.val = false;
end

% degeneracy
if (~isempty(P1.fullDim.val) && ~P1.fullDim.val) ...
        || (~isempty(P2.fullDim.val) && ~P2.fullDim.val)
    P_out.fullDim.val = false;
elseif (~isempty(P1.fullDim.val) && P1.fullDim.val) ...
        && (~isempty(P2.fullDim.val) && P2.fullDim.val)
    P_out.fullDim.val = true;
end

% minimal representations
if (~isempty(P1.minHRep.val) && P1.minHRep.val) ...
        || (~isempty(P2.minHRep.val) && P2.minHRep.val)
    P_out.minHRep.val = true;
elseif (~isempty(P1.minHRep.val) && ~P1.minHRep.val) ...
        && (~isempty(P2.minHRep.val) && ~P2.minHRep.val)
    P_out.minHRep.val = false;
end
if (~isempty(P1.minVRep.val) && P1.minVRep.val) ...
        || (~isempty(P2.minVRep.val) && P2.minVRep.val)
    P_out.minVRep.val = true;
elseif (~isempty(P1.minVRep.val) && ~P1.minVRep.val) ...
        && (~isempty(P2.minVRep.val) && ~P2.minVRep.val)
    P_out.minVRep.val = false;
end

end

function P_out = aux_cartProd_polytope_point(P,p,order)
% Cartesian product for polytope and point
% helper variable: order = true: P x p; order = false: p x P

% read out dimension
n = size(p,1);
% init equality constraints
Ae = eye(n);
be = p;

% instantiate polytope
if order
    % polytope x point
    P_out = cartProd_(P,polytope(zeros(0,n),[],Ae,be),'exact');
else
    % point x polytope
    P_out = cartProd_(polytope(zeros(0,n),[],Ae,be),p,'exact');
end
% resultng polytope is degenerate
P_out.fullDim.val = false;
% set is (non)empty if the polytope is (non)empty
P_out.emptySet.val = P.emptySet.val;
% set is (un)bounded if the polytope is (un)bounded
P_out.bounded.val = P.bounded.val;

end

% ------------------------------ END OF CODE ------------------------------
