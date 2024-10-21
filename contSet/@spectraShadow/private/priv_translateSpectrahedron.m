function SpS_out = priv_translateSpectrahedron(SpS,t)
% priv_translateSpectrahedron - translates a 'pure' spectrahedron by
%    modifying its coefficient matrices A0,A1,...,An, not its center vector.
%    Important: This assumes that SpS is a 'pure' spectrahedron, i.e.,
%    SpS.G = eye(n), SpS.c = zeros([n 1]) (this is not checked!)
%
% Syntax:
%    SpS = priv_translateSpectrahedron(SpS,t)
%
% Inputs:
%    SpS - spectraShadow object representing a 'pure' spectrahedron
%    t - numerical vector
%
% Outputs:
%    SpS - translated spectraShadow object
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Adrian Kulmburg
% Written:       27-August-2023 
% Last update:   ---    
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

[A0,Ai] = priv_getCoeffMatrices(SpS);
for i=1:dim(SpS)
    A0 = A0 - t(i) * Ai{i};
end
SpS_out = spectraShadow([A0 cat(2, Ai{:})]);

% properties
SpS_out.bounded.val = SpS.bounded.val;
SpS_out.emptySet.val = SpS.emptySet.val;
SpS_out.fullDim.val = SpS.fullDim.val;
if ~isempty(SpS.center.val)
    SpS_out.center.val = SpS_out.center.val + t;
end

% ------------------------------ END OF CODE ------------------------------
