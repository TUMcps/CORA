function res = testLong_polyZonotope_jacobian
% testLong_polyZonotope_jacobian - unit test function for computing
%    the jacobian (complete derivative) of a polynomial zonotope
%
% Syntax:
%    res = testLong_polyZonotope_jacobian
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann
% Written:       24-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
nTests = 10;

for i=1:nTests
    %% analytic test
    n = randi(30);
    ng = randi([3,10]);
    nf = ng - 1;
    pZ = noIndep(polyZonotope.generateRandom('Dimension',n,...
        'NrGenerators',ng,'NrFactors',nf));
    x = sym('x',[size(pZ.E,1),1],'real');
    ne = length(pZ.id);
    ind_diff = ismember(pZ.id,unique(randi(ne-1,ne-1,1)));
    
    % compute "analytic" result
    f = fhandle(pZ);
    f_sym = f(x);
    jac_sym = jacobian(f_sym,x(ind_diff));
    
    % compute jacobian of pZ
    pZ_diff_c = jacobian(pZ,pZ.id(ind_diff));
    for j=1:length(pZ_diff_c)
        fj = fhandle(restoreId(pZ_diff_c{j},pZ.id),{pZ.id});
        valj = fj(x)-jac_sym(:,j);
        for k=1:n
            ck = coeffs(valj(k),x);
            if any(ck>1e-8)
                res = false;
                return;
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
