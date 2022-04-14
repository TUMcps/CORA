function res = testLongDuration_polyZonotope_fhandle
% testLongDuration_polyZonotope_fhandle - unit test function for computing
%    the function handle of a polynomial zonotope
%
% Syntax:  
%    res = testLongDuration_polyZonotope_fhandle
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
nTests = 50;
for i=1:nTests
    n = randi(30);
    pZ = noIndep(polyZonotope.generateRandom(n));
    ne = length(pZ.id);
    f = fhandle(pZ);
    N = 2*n;
    B = 2*rand(ne,N)-1;
    for j=1:N
        fval_res = resolve(pZ,B(:,j));
        fval = f(B(:,j));
        f_abs = max(abs(fval_res),abs(fval));
        f_abs(f_abs==0) = 1;
        if any(abs(fval_res-fval)./f_abs>1e-6)
            res = false;
            break;
        end
    end
    if ~res
        break;
    end
end

if ~res
    disp('testLongDuration_polyZonotope_fhandle failed');
else
    disp('testLongDuration_polyZonotope_fhandle successful');
end
%------------- END OF CODE --------------