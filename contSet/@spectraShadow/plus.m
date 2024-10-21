function S_out = plus(SpS,S)
% plus - Overloaded '+' operator for the minkowski addition of a
%    spectrahedral shadow and a vector, or a spectrahedral shadow and a
%    contSet object
%
% Syntax:
%    S_out = SpS + S
%    S_out = plus(SpS,S)
%
% Inputs:
%    SpS - spectraShadow object, numeric
%    S - contSet object, numeric vector
%
% Outputs:
%    S_out - spectraShadow object, result of the addition of SpS and obj
%
% Example:
%    SpS = spectraShadow([1 0 1 0;0 1 0 -1]);
%    SpS_times_2 = SpS + SpS;
%    SpS_translated = SpS + 5;
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Maximilian Perschl, Adrian Kulmburg
% Written:       24-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% ensure that numeric is second input argument
[SpS,S] = reorderNumeric(SpS,S);

% check dimensions
equalDimCheck(SpS,S);

% addition between spectrahedral shadow and numerical vector
if isnumeric(S) && iscolumn(S)
    S_out = aux_plus_vector(SpS,S);
    return
end

% minkowski addition between two spectrahedra
if isa(S,'spectraShadow')
    S_out = aux_plus_spectraShadow(SpS,S);
    return
end

if isa(S,'contSet')
    try
        S = spectraShadow(S);
    catch ME
        throw(CORAerror('CORA:noExactAlg',SpS,S));
    end
    S_out = aux_plus_spectraShadow(SpS,S);
    return
end

throw(CORAerror('CORA:noops',SpS,S));

end


% Auxiliary functions -----------------------------------------------------

function S_out = aux_plus_spectraShadow(SpS,S)

c_1 = SpS.c;
c_2 = S.c;

G_1 = SpS.G;
G_2 = S.G;

[A0_1,Ai_1] = priv_getCoeffMatrices(SpS);
[A0_2,Ai_2] = priv_getCoeffMatrices(S);

m_1 = size(G_1,2);
m_2 = size(G_2,2);
k_1 = size(A0_1,1);
k_2 = size(A0_2,1);

c = c_1 + c_2;
G = [G_1 G_2];
A0 = blkdiag(A0_1,A0_2);
Ai = cell([1 m_1+m_2]);
for i=1:m_1
    Ai{i} = blkdiag(Ai_1{i},sparse(k_2,k_2));
end
for i=1:m_2
    Ai{m_1+i} = blkdiag(sparse(k_1,k_1),Ai_2{i});
end

S_out = spectraShadow([A0 cat(2, Ai{:})],c,G);

% Adding some extra information
if ~isempty(SpS.bounded.val) && ~isempty(S.bounded.val)
    if SpS.bounded.val && S.bounded.val
        S_out.bounded.val = true;
    else
        S_out.bounded.val = false;
    end
end
if (~isempty(SpS.fullDim.val) && SpS.fullDim.val) ...
        || (~isempty(S.fullDim.val) && S.fullDim.val)
    S_out.fullDim.val = true;
end
if (~isempty(SpS.emptySet.val) && SpS.emptySet.val) ...
        || (~isempty(S.emptySet.val) && S.emptySet.val)
    S_out.emptySet.val = true;
elseif (~isempty(SpS.emptySet.val) && ~SpS.emptySet.val) ...
        && (~isempty(S.emptySet.val) && ~S.emptySet.val)
    S_out.emptySet.val = false;
end

if ~isempty(SpS.center.val) && ~isempty(S.center.val)
    S_out.center.val = SpS.center.val + S.center.val;
end

end

function S_out = aux_plus_vector(SpS,S)

% get parameters of S
A = SpS.A; G = SpS.G; c = SpS.c;

% compute parameters of new spectrahedron
if isemptyobject(SpS)
    c_new = [];
else
    c_new = S+c;
end

% construct spectrahedral shadow
S_out = spectraShadow(A,c_new,G);

S_out.bounded.val = SpS.bounded.val;
S_out.fullDim.val = SpS.fullDim.val;
S_out.emptySet.val = SpS.emptySet.val;

if ~isempty(SpS.center.val)
    S_out.center.val = SpS.center.val + S;
end

end

% ------------------------------ END OF CODE ------------------------------
