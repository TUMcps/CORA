function SpS_out = cartProd_(SpS,S,type,varargin)
% cartProd_ - computes the Cartesian product of a spectrahedral shadow and
%    a contSet object
%
% Syntax:
%    SpS_out = cartProd_(SpS,S)
%
% Inputs:
%    SpS - spectraShadow object
%    S - contSet object
%    type - 'exact' (other methods not yet implemented)
%
% Outputs:
%    SpS_out - spectraShadow object
%
% Example: 
%    SpS_segment = spectraShadow([1 0 1 0;0 1 0 -1]);
%    SpS_box = cartProd(SpS_segment, SpS_segment);
%    plot(SpS_box)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/cartProd

% Authors:       Maximilian Perschl, Adrian Kulmburg
% Written:       12-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if ~isa(S,'spectraShadow')
    % Try to recast the set as spectrahedral shadow
    try
        S = spectraShadow(S);
    catch ME
        throw(CORAerror("CORA:noops",SpS,S));
    end
end

% dimensions of both sets
d = dim(SpS);
e = dim(S);

% If one of the sets is known to represent the empty set, we can get it 
% done a bit faster
if (~isempty(SpS.emptySet.val) && SpS.emptySet.val) ...
        || (~isempty(S.emptySet.val) && S.emptySet.val)
    A = zeros([1 d+e+1]);
    A(1) = -1;
    SpS_out = spectraShadow(A);
    SpS_out.emptySet.val = true;
    SpS_out.bounded.val = true;
    SpS_out.fullDim.val = false;
    SpS_out.center = [];
    return
end

% cartesian product of S1 and R^e
SpS_e = aux_rightCartProdSpace(SpS,e);

% cartesian product of S2 and R^d
S_d = aux_leftCartProdSpace(S,d);

% result is the intersection of the previous results
SpS_out = SpS_e & S_d;

% If both sets are known to be bounded, we can deduce the same about the
% cartProd
if (~isempty(SpS.bounded.val) && SpS.bounded.val) ...
        && (~isempty(S.bounded.val) && S.bounded.val)
    SpS_out.bounded.val = true;
end
% Same for full-dimensionality
if (~isempty(SpS.fullDim.val) && SpS.fullDim.val) ...
        && (~isempty(S.fullDim.val) && S.fullDim.val)
    SpS_out.fullDim.val = true;
end
% And similarly for the center
if ~isempty(SpS.center.val) ...
        && ~isempty(S.center.val)
    SpS_out.center.val = [SpS.center.val;S.center.val];
end

end


% Auxiliary functions -----------------------------------------------------

function SpS_Rd = aux_rightCartProdSpace(SpS,d)
% compute the cartesian product S x R^d
    % create new projection
    G = SpS.G;
    c = SpS.c;
    
    new_G = [G sparse(dim(SpS),d); sparse(d,size(G,2)) speye(d)];
    new_c = [c;zeros(d,1)];

    % create zero-matrix coefficients and concatenate them to the previous
    % coefficients
    A = SpS.A;
    k = size(A,1);
    new_A = [A sparse(k,d*k)];

    % new spectrahedral shadow
    SpS_Rd = spectraShadow(new_A,new_c,new_G);
end

function SpS_Rd = aux_leftCartProdSpace(SpS,d)
% compute the cartesian product R^d x S
    % create new projection
    G = SpS.G;
    c = SpS.c;
    
    new_G = [speye(d) sparse(d,size(G,2)); sparse(size(G,1),d) G];
    new_c = [zeros(d,1);c];

    % create zero-matrix coefficients and concatenate them to the previous
    % coefficients
    A = SpS.A;
    k = size(A,1);
    new_A = [A(:,1:k) sparse(k,d*k) A(:,k+1:end)];

    % new spectrahedral shadow
    SpS_Rd = spectraShadow(new_A,new_c,new_G);
end

% ------------------------------ END OF CODE ------------------------------
