function P_out = box(P)
% box - computes an enclosing axis-aligned box represented as a polytope in
%    halfspace representation
%
% Syntax:
%    P_out = box(P)
%
% Inputs:
%    P - polytope object
%
% Outputs:
%    P_out - polytope object 
%
% Example:
%    A = [1 2; -2 1; -2 -2; 3 -1];
%    b = ones(4,1);
%    P = polytope(A,b);
%    B = box(P);
%
%    figure; hold on;
%    plot(P);
%    plot(B,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Viktor Kotsev, Mark Wetzlinger
% Written:       16-May-2022
% Last update:   14-December-2022 (MW, unbounded case, MOSEK support, remove equality constraints)
% Last revision: 12-July-2024 (MW, refactor)

% ------------------------------ BEGIN CODE -------------------------------

% fullspace case
if representsa_(P,'fullspace',0)
    P_out = polytope.Inf(dim(P)); return
end

% the computation of the box outer approximation is much faster for the
% vertex representation, so we first check for that
if P.isVRep.val
    % --- V representation: take the min and max along each dimension
    P_out = aux_box_V(P);
else
    % --- H representation: compute the support function in the direction
    % of all 2n plus/minus axis-aligned basis vectors
    P_out = aux_box_H(P);
end

end


% Auxiliary functions -----------------------------------------------------

function P_out = aux_box_V(P)
% computation of box enclosure for a polytope in vertex representation

% compute halfspace representation
[A,b,empty,fullDim,bounded] = priv_box_V(P.V_.val,dim(P));
if empty
    % add properties to input polytope
    P.emptySet.val = true;
    P.bounded.val = true;
    P.fullDim.val = false;
    P.V_.val = zeros(dim(P),0);
    P.isVRep.val = true;
    P.minVRep.val = true;
    % output empty polytope
    P_out = polytope.empty(dim(P));
    return
end
% instantiate polytope (note that this eliminates all constraints where the
% offset b is +-Inf)
P_out = polytope(A,b);

% set properties
P_out.minHRep.val = true;
P_out.emptySet.val = empty;
P_out.fullDim.val = fullDim;
P_out.bounded.val = bounded;

end

function P_out = aux_box_H(P)
% computation of box enclosure for a polytope in halfspace representation

% compute halfspace representation
[A,b,empty,fullDim,bounded] = priv_box_H(P.A_.val,P.b_.val,...
    P.Ae_.val,P.be_.val,dim(P));

if empty
    % add properties to input polytope
    P.emptySet.val = true;
    P.bounded.val = true;
    P.fullDim.val = false;
    P.V_.val = zeros(dim(P),0);
    P.isVRep.val = true;
    P.minVRep.val = true;
    % output empty polytope
    P_out = polytope.empty(dim(P));
    return
end

% set properties: input polytope
P.emptySet.val = false;
P.fullDim = fullDim;
P.bounded.val = bounded;

% construct box
P_out = polytope(A,b);
% set properties: output polytope
P_out.minHRep.val = true;
P_out.emptySet.val = empty;
P_out.fullDim.val = fullDim;
P_out.bounded.val = bounded;

end

% ------------------------------ END OF CODE ------------------------------
