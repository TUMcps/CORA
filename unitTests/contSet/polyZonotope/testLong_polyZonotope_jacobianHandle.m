function res = testLong_polyZonotope_jacobianHandle
% testLong_polyZonotope_jacobianHandle - unit test function for
%    computing the fhandle of the jacobian of a polyZonotope
%
% Syntax:
%    res = testLong_polyZonotope_jacobianHandle
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
    
    % compute jacobian handle of pZ
    f_jac = jacobianHandle(pZ,pZ.id(ind_diff),pZ.id(~ind_diff));
    F_val = simplify(jac_sym-f_jac(x(ind_diff),x(~ind_diff)));
    f_val = F_val(:);
    for j=1:length(f_val)
        cj = coeffs(f_val(j),x);
        if any(cj>1e-8)
            res = false;
            return;
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
