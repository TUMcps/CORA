function res = priv_mtimesIsomorphism(SpS,M)
% priv_mtimesIsomorphism - multiplication of a matrix M with a
%    spectrahedron S, but only in the case where M is an isomorphism, i.e.,
%    is bijective (this is not checked, since it's a private function).
%    Such a transformation can leave the generator matrix G invariant, and
%    only modify the coefficient matrices in A instead.
%    IMPORTANT: This function will ALSO modify the existential sum
%    representation.
%
% Syntax:
%    res = priv_mtimesIsomorphism(SpS,M)
%
% Inputs:
%    SpS - spectraShadow object
%    M - Quadratic, full rank (!) matrix
%
% Outputs:
%    res - the resulting spectrahedron
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Adrian Kulmburg
% Written:       02-August-2023 
% Last update:   ---    
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

c = SpS.c;
G = SpS.G;
A = SpS.A;
ESumRep = SpS.ESumRep.val;

n = dim(SpS);
m = size(G,2);

Minv = inv(M);

% We now compute the new existential sum representation
A_ESumRep = ESumRep{1};
B_ESumRep = ESumRep{2};

if isempty(A_ESumRep)
    A_ESumRep_new = [];
else
    % Dirty cast to quickly get the coeff matrices of A
    [A0_ESumRep, Ai_ESumRep] = priv_getCoeffMatrices(spectraShadow(A_ESumRep));
    % Now, transform them as before
    A0_ESumRep_new = A0_ESumRep;
    Ai_ESumRep_new = cell([1 n]);
    for i=1:n
        Ai_ESumRep_new{i} = 0;
        for j=1:n
            % Manually multiply the coefficient matrices
            Ai_ESumRep_new{i} = Ai_ESumRep_new{i} + Minv(j,i) * Ai_ESumRep{j};
        end
    end
    A_ESumRep_new = [A0_ESumRep_new cat(2,Ai_ESumRep_new{:})];
end

% Since we have an explicit formula for the ESumRep, we can save it
ESumRep_new = {A_ESumRep_new B_ESumRep};

res = spectraShadow(A, M*c, M*G);
res.ESumRep.val = ESumRep_new;

% Set properties
res.bounded.val = SpS.bounded.val;
res.fullDim.val = SpS.fullDim.val;
res.emptySet.val = SpS.emptySet.val;
if ~isempty(SpS.center.val)
    res.center.val = M * SpS.center.val;
end

% ------------------------------ END OF CODE ------------------------------
