function res = taylm(obj)
% taylm - enclose a polyZonotope object with a Talyor model
%
% Syntax:  
%    res = taylm(obj)
%
% Inputs:
%    obj - polyZonotope object
%
% Outputs:
%    res - taylm object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, zonotope

% Author:       Niklas Kochdumper
% Written:      13-August-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
    % introduce independent factors as new dependent factors
    G = [obj.G,obj.Grest];
    expMat = blkdiag(obj.expMat,eye(size(obj.Grest,2)));

    % create Taylor model for factors
    p = length(obj.id) + size(obj.Grest,2);
    int = interval(-ones(p,1),ones(p,1));
    
    tay = taylm(int);
    
    % convert polyZonotope object to taylor model
    res = obj.c;
    
    for i = 1:size(G,2)
        temp = 1;
        for j = 1:size(expMat,1)
           temp = temp * tay(j)^expMat(j,i); 
        end
        res = res + G(:,i) * temp;
    end

%------------- END OF CODE --------------