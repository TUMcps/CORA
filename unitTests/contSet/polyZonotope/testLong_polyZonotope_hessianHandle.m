function res = testLong_polyZonotope_hessianHandle
% testLong_polyZonotope_hessianHandle - unit test function for
%    computing the fhandle of the hessian of a 1D polyZonotope
%
% Syntax:
%    res = testLong_polyZonotope_hessianHandle
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

    % random dimension
    n = randi(1);

    % random number of generators
    ng = randi([3,10]);

    % number of factors
    nf = ng - 1;

    % instantiate polynomial zonotope
    pZ = noIndep(polyZonotope.generateRandom('Dimension',n,...
        'NrGenerators',ng,'NrFactors',nf));
    
    x = sym('x',[size(pZ.E,1),1],'real');
    ne = length(pZ.id);
    ind_diff = ismember(pZ.id,unique(randi(ne-1,ne-1,1)));
    
    % compute "analytic" result
    f = fhandle(pZ);
    f_sym = f(x);
    hess_sym = hessian(f_sym,x(ind_diff));
    
    % compute hessian handle of pZ
    [f_hess,H_str] = hessianHandle(pZ,pZ.id(ind_diff),pZ.id(~ind_diff));
    
    % check structure
    if ~isempty(nonzeros(hess_sym(~H_str)))
        res = false;
        return;
    end
    
    % check function handle
    F_val = simplify(hess_sym-f_hess(x(ind_diff),x(~ind_diff)));
    if isempty(nonzeros(F_val))
        continue;
    end
    f_val = F_val(:);
    for j=1:length(f_val)
        cj = abs(double(coeffs(f_val(j),x)));
        if any(cj > 0 & ~withinTol(cj,0,1e-8))
            res = false;
            return;
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
