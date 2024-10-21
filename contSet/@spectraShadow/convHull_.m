function S_out = convHull_(SpS,S,varargin)
% convHull_ - returns the convex hull over the union of a spectrahedral
%    shadow and a contSet
%
% Syntax:
%    SpS_out = convHull_(SpS,S)
%
% Inputs:
%    SpS - spectraShadow object
%    S - contSet object
%
% Outputs:
%    S_out - convex hull over the union of SpS and S
%
% Example:
%    A0 = eye(4);
%    A1 = [-1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0];
%    A2 = [0 0 0 0;0 0 0 0;0 0 -1 0;0 0 0 1];
%    SpS = spectraShadow([A0 A1 A2]);
%    SpS_convHull = convHull(SpS,SpS+[1;1]);
%    
%    hold on
%    plot(SpS)
%    plot(SpS+[1;1])
%    plot(SpS_convHull)
% 
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/convHull

% Authors:       Maximilian Perschl, Adrian Kulmburg
% Written:       16-May-2023 
% Last update:   ---
% Last revision: 29-September-2024 (MW, integrate precedence)

% ------------------------------ BEGIN CODE -------------------------------

% spectraShadow is already convex
if nargin == 1
    S_out = SpS;
    return;
end

% ensure that numeric is second input argument
[SpS,S] = reorderNumeric(SpS,S);

% check dimensions
equalDimCheck(SpS,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < SpS.precedence
    S_out = convHull(S,SpS,varargin{:});
    return
end

% all cases: convert to spectraShadow
try
    S = spectraShadow(S);
catch ME
    throw(CORAerror('CORA:noops',SpS,S));
end

% empty cases
if isemptyobject(SpS) 
    S_out = S;
    return;
end

if isemptyobject(S)
    S_out = SpS;
    return;
end

% construct R^n x {1}
n = dim(SpS);
c_one = [zeros(n,1);1];
Rn_x_1 = spectraShadow(sparse(1,n+1),c_one,[speye(n);sparse(1,n)]);

% Construct singleton, i.e., {1}
singleton = spectraShadow(sparse([1 0 -1 0;0 -1 0 1]));
% compute S1 x {1} and S2 x {1}
S1_x_1 = cartProd_(SpS,singleton,'exact');
S2_x_2 = cartProd_(S,singleton,'exact');

% compute conic hulls
cone_S1 = aux_conicHull(S1_x_1);
cone_S2 = aux_conicHull(S2_x_2);

cone_sum = cone_S1 + cone_S2;

% compute final result
S_out = [speye(n) sparse(n, 1)] * (cone_sum & Rn_x_1);

% Filling out additional information
% S1 or S2 full Dim -> S full Dim
if (~isempty(SpS.fullDim.val) && SpS.fullDim.val) || (~isempty(S.fullDim.val) && S.fullDim.val)
    S_out.fullDim.val = true;
end

if ~isempty(SpS.bounded.val) && ~isempty(S.bounded.val)
    if SpS.bounded.val && S.bounded.val
        % S1 and S2 bounded -> S bounded
        S_out.bounded.val = true;
    else
        S_out.bounded.val = false;
    end
end

if ~isempty(SpS.center.val) && ~isempty(S.center.val)
    % New center is average of the centers of each set
    S_out.center.val = 0.5 * (SpS.center.val + S.center.val);
end

end


% Auxiliary functions -----------------------------------------------------

function coneHull = aux_conicHull(SpS)
% Compute the conic hull of the spectrahedral shadow SpS
% For a spectrahedral shadow given in existential sum representation
% SpS = {x|\exists y, B0 + x1*B1 + ... + xn*Bn + y1*C1 + ... + yl*Cl >= 0}
% the conic hull is given as
% cone(SpS) = {x | \exists y, \exists s,t
% s*B0+x1*B1+...+xn*Bn+y1*C1+...+yl*Cl>=0 and
% xi*[0 1;1 0]+s*[1 0;0 0]+t*[0 0;0 1]>=0 for all i}

generateESumRep(SpS);
ESumRep = SpS.ESumRep.val;
B = ESumRep{1};
% We separate B
[B0, Bi] = priv_getCoeffMatrices(spectraShadow(B));
C = ESumRep{2};

% Relevant sizes
n = dim(SpS);
k = size(B0,1);

% First, we begin by constructing the new B.
% In total, we need to add 2*n new constraints, so our new B0 is
B0_new = sparse(k+2*n,k+2*n);
% We continue with the new Bi; for each individual one, we have to
% block-diagonalize one matrix [0 1;1 0], and n-1 2x2-zero matrices:
Bi_new = cell([1 n]);
for i=1:n
    Bi_new{i} = blkdiag(Bi{i}, sparse(2*(i-1),2*(i-1)),sparse([0 1;1 0]),sparse(2*(n-i),2*(n-i)));
end

B_new = [B0_new cat(2, Bi_new{:})];

% We now construct the new C. The first coefficient matrices, corresponding
% to the yi, just need to be expanded vertically:
L = size(C,2)/k;
C_start = cell([1 L]);
for i=1:L
    C_start{i} = blkdiag(C(:,(i-1)*k+1:i*k), sparse(2*n, 2*n));
end
C_new = cat(2, C_start{:});
% Next up are the coefficient matrices for the s:
C_s = blkdiag(B0, kron(speye(n), sparse([1 0;0 0])));
C_new = [C_new C_s];
% And finally, we build the coefficient matrices for the t:
C_t = blkdiag(sparse(k,k), kron(speye(n), sparse([0 0;0 1])));
C_new = [C_new C_t];

ESumRep_new = {B_new C_new};

coneHull = spectraShadow(ESumRep_new);

end

% ------------------------------ END OF CODE ------------------------------
