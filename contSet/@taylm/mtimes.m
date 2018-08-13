function res = mtimes(factor1, factor2)
% mtimes - Overloaded '*' operator for a Taylor model
%
% Syntax:  
%    res = mtimes(factor1, factor2)
%
% Inputs:
%    factor1 and factor2 - taylm objects
%
% Outputs:
%    res - a taylm object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm

% Author:       Dmitry Grebenyuk
% Written:      20-August-2017
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE -------------
[m1,n1] = size(factor1);
[m2,n2] = size(factor2);
if n1 == m2
    res = taylm();
    res(m1,n2) = taylm();
    T = 0;

    for i = 1:m1
        for j = 1:n2
            for k = 1:n1
                T = T + factor1(i,k) .* factor2(k,j);
            end
            res(i,j) = T;
            T = 0;
        end
    end
    
elseif (isa(factor1, 'double') || isa(factor2, 'double')) && (isscalar(factor1) || isscalar(factor2))
    res = factor1 .* factor2;
else
    error('Matrix dimensions do not agree!');
end
%------------- END OF CODE --------------
