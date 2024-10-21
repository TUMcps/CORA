function SpS_out = and_(SpS,S,varargin)
% and_ - Overloaded '&' operator for the intersection of a spectrahedral
%    shadow and a contSet object
%
% Syntax:
%    SpS_out = and_(SpS,S)
%    SpS_out = and_(SpS,S)
%
% Inputs:
%    SpS - spectraShadow object
%    S - contSet object
%
% Outputs:
%    SpS_out - spectraShadow object
%
% Example:
%    A0 = eye(4);
%    A1 = [-1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0];
%    A2 = [0 0 0 0;0 0 0 0;0 0 -1 0;0 0 0 1];
%    SpS = spectraShadow([A0 A1 A2]);
%    SpS_rectangle = SpS & (SpS + [0.5;0]);
%    hold on
%    plot(SpS)
%    plot(SpS+[0.5;0])
%    plot(SpS_rectangle)
%    
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/and

% Authors:       Maximilian Perschl, Adrian Kulmburg
% Written:       24-April-2023
% Last update:   ---
% Last revision: 28-September-2024 (MW, integrate precedence)

% ------------------------------ BEGIN CODE -------------------------------

% call function with lower precedence
if isa(S,'contSet') && S.precedence < SpS.precedence
    SpS_out = and_(S,SpS,varargin{:});
    return
end

% spectraShadow
if isa(S,'spectraShadow')
    SpS_out = aux_and_spectraShadow(SpS,S);
    return
end

% all other cases: conversion to spectraShadow
try
    S = spectraShadow(S);
catch ME
    throw(CORAerror('CORA:noops',SpS,S));
end
SpS_out = aux_and_spectraShadow(SpS,S);

end


% Auxiliary functions -----------------------------------------------------

function SpS_out = aux_and_spectraShadow(SpS,S)

% If both sets are known to represent the empty set, we can get it done a
% bit faster
if ~isempty(SpS.emptySet.val) && ~isempty(S.emptySet.val)
    if SpS.emptySet.val && S.emptySet.val
        A = zeros([1 dim(SpS)+1]);
        A(1) = -1;
        SpS_out = spectraShadow(A);
        SpS_out.emptySet.val = true;
        SpS_out.bounded.val = true;
        SpS_out.fullDim.val = false;
        return
    end
end


% generate suitable representation of both spectrahedra
generateESumRep(SpS);
generateESumRep(S);

% We extract the existential sum representation
B1 = SpS.ESumRep.val{1};
C1 = SpS.ESumRep.val{2};
B2 = S.ESumRep.val{1};
C2 = S.ESumRep.val{2};

k1 = size(B1,1);
k2 = size(B2,1);

% We need to extract the coefficient matrices from B1 and B2; this can most
% easily be done by defining a spectraShadow with B1 and B2 respectively,
% and extracting the coefficient matrices:
[constant_B1, dependent_B1] = priv_getCoeffMatrices(spectraShadow(B1));
[constant_B2, dependent_B2] = priv_getCoeffMatrices(spectraShadow(B2));

% For the new ESumRep{1}, it suffices to take the block-diagonal
% combination of the matrices for B1 and B2:
constant_new = blkdiag(constant_B1, constant_B2);
dependent_new = cell([1 size(dependent_B1,2)]);
for i=1:size(dependent_B1,2)
    dependent_new{i} = blkdiag(dependent_B1{i}, dependent_B2{i});
end

ESumRep_new{1} = [constant_new cat(2,dependent_new{:})];

% For the new ESumRep{2}, it is quite similar, except that we need to
% block-diagonalize with zero-matrices
% Define some zero-matrices
zeros1 = sparse(k1,k1);
zeros2 = sparse(k2,k2);

% Some helpful lengths
m1 = size(C1,2)/k1;
m2 = size(C2,2)/k2;

C_new = cell([1 m1+m2]);
for i=1:m1
    C1_part = C1(:,k1*(i-1)+1:k1*i);
    C_new{i} = blkdiag(C1_part, zeros2);
end
for i=1:m2
    C2_part = C2(:,k2*(i-1)+1:k2*i);
    C_new{m1+i} = blkdiag(zeros1, C2_part);
end
C_new = cat(2, C_new{:});
ESumRep_new{2} = C_new;

% generate spectraShadow
SpS_out = spectraShadow(ESumRep_new);

% If both sets are known to be bounded, we can deduce the same for S
if ~isempty(SpS.bounded.val) && ~isempty(S.bounded.val)
    if SpS.bounded.val && S.bounded.val
        SpS_out.bounded.val = true;
    end
end

end

% ------------------------------ END OF CODE ------------------------------
