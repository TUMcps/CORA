function res = testLongDuration_polyZonotope_hessianHandle
% testLongDuration_polyZonotope_hessianHandle - unit test function for
%    computing the fhandle of the hessian of a 1D polyZonotope
%
% Syntax:  
%    res = testLongDuration_polyZonotope_hessianHandle
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      24-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;
nTests = 10;

for i=1:nTests
    %% analytic test
    n = randi(1);
    ng = randi([3,10]);
    nf = randi([2,10]);
    pZ = noIndep(polyZonotope.generateRandom(n,ng,nf));
    x = sym('x',[size(pZ.expMat,1),1],'real');
    ne = length(pZ.id);
    ind_diff = ismember(pZ.id,unique(randi(ne-1,ne-1,1)));
    
    % compute "analytic" result
    f = fhandle(pZ);
    f_sym = f(x);
    hess_sym = hessian(f_sym,x(ind_diff));
    
    % compute hessian handle of pZ
    [f_hess,H_str] = hessianHandle(pZ,pZ.id(ind_diff));
    
    % check structure
    if ~isempty(nonzeros(hess_sym(~H_str)))
        res = false;
        break;
    end
    
    % check function handle
    F_val = simplify(hess_sym-f_hess(x(ind_diff),x(~ind_diff)));
    if isempty(nonzeros(F_val))
        continue;
    end
    f_val = F_val(:);
    for j=1:length(f_val)
        cj = abs(double(coeffs(f_val(j),x)));
        if any(cj>1e-8)
            res = false;
            break;
        end
    end
    if ~res
        break;
    end
end
if ~res
    disp('testLongDuration_polyZonotope_hessianHandle failed');
else
    disp('testLongDuration_polyZonotope_hessianHandle successful');
end
%------------- END OF CODE --------------