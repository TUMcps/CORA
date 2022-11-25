function res = testLongDuration_polyZonotope_resolve
% testLongDuration_polyZonotope_resolve - unit test function for
%    (partially) resolving a polynomial zonotope
%
% Syntax:  
%    res = testLongDuration_polyZonotope_resolve
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
    n = randi(30);
    ng = randi([3,10]);
    nf = randi([2,10]);
    pZ = noIndep(polyZonotope.generateRandom(n,ng,nf));
    x = sym('x',[size(pZ.expMat,1),1],'real');
    ne = length(pZ.id);
    ind_res = ismember(pZ.id,unique(randi(ne-1,ne-1,1)));
    val_res = 2*rand(sum(ind_res),1)-1;
    f = fhandle(pZ,{pZ.id(~ind_res),pZ.id(ind_res)});
    f_val = f(x(~ind_res),val_res);

    % resolve
    pZ_res = resolve(pZ,val_res,pZ.id(ind_res));
    f_res = fhandle(pZ_res);
    fres_val = f_res(x(~ind_res));
    sym_res = simplify(f_val-fres_val);
    for j=1:n
        cj = coeffs(sym_res(j),x(~ind_res));
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
    disp('testLongDuration_polyZonotope_resolve failed');
else
    disp('testLongDuration_polyZonotope_resolve successful');
end
%------------- END OF CODE --------------