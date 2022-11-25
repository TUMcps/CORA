function res = isZonotope(pZ)
% isZonotope - Checks if a polynomial zonotope represents a zonotope
%
% Syntax:  
%    res = isZonotope(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    res - 1 if set is a zonotope, 0 if not
%
% Example: 
%    pZ1 = polyZonotope([-0.5;0],[1.5 -0.5;-0.5 -2],[],[0 1;1 0]);
%    pZ2 = polyZonotope([-0.5;0],[-0.5 -0.5;0.5 -2],[],[1 1;0 1]);
%   
%    isZonotope(pZ1)
%    isZonotope(pZ2)
%
%    figure
%    hold on
%    plot(pZ1,[1,2],'b','Filled',true,'EdgeColor','none','Splits',10);
%
%    figure
%    hold on
%    plot(pZ2,[1,2],'r','Filled',true,'EdgeColor','none','Splits',10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope, isPolytope

% Author:       Niklas Kochdumper
% Written:      14-August-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    res = 0;

    % remove redundant exponent vectors
    [expMat,G] = removeRedundantExponents(pZ.expMat,pZ.G);

    % check matrix dimensions
    if size(expMat,1) ~= size(G,2)
        return;
    end

    % sort exponent matrix rows
    expMat = sortrows(expMat,'descend');

    % check if exponent matrix is the identity matrix
    if sum(sum(abs(expMat-diag(diag(expMat))))) == 0
        res = 1;
    end
end

%------------- END OF CODE --------------
