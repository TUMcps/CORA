function pZ = polyZonotope(cPZ)
% polyZonotope - enclose constrained polynomial zonotope with a polynomial
%                zonotope
%
% Syntax:  
%    pZ = polyZonotope(cPZ)
%
% Inputs:
%    cPZ - conPolyZono object 
%
% Outputs:
%    pZ - polyZonotope object enclosing cPZ
%
% Example:  
%    c = [0;0];
%    G = [1 0 1 -1; 0 2 1 2];
%    expMat = [1 2 1 0; 0 0 1 2; 0 0 0 0];
%    A = [1 1 0.5];
%    b = 0.5;
%    expMat_ = [0 1 0;1 0 0; 0 0 1];
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_);
%
%    pZ = polyZonotope(cPZ);
%   
%    figure; hold on;
%    plot(cPZ,[1,2],'r','Splits',20,'Filled',true,'EdgeColor','none');
%    plot(pZ);
%
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope, interval, zonotope

% Author:       Niklas Kochdumper
% Written:      04-February-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % remove all constraints
    cPZ = reduceConstraints(cPZ,0);
    
    % represent the resulting object as a polynomial zonotope
    pZ = polyZonotope(cPZ.c,cPZ.G,cPZ.Grest,cPZ.expMat,cPZ.id);
    
end
    
%------------- END OF CODE --------------