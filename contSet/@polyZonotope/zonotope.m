function Z = zonotope(pZ)
% zonotope - computes an enclosing zonotope of the polynomial zonotope
%
% Syntax:  
%    Z = zonotope(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    pZ = polyZonotope([0;0],[2 0 1 1;0 2 1 2],[0;0],[1 0 3 1;0 1 0 2]);
%    zono = zonotope(pZ);
%
%    figure
%    hold on
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%    plot(zono,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, mptPolytope

% Author:       Niklas Kochdumper
% Written:      24-March-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

if ~isempty(pZ.G)
    
    % determine dependent generators with exponents that are all even
    temp = prod(ones(size(pZ.expMat))-mod(pZ.expMat,2),1);
    Gquad = pZ.G(:,temp == 1);

    % compute zonotope parameter
    c = pZ.c + 0.5 * sum(Gquad,2);
    G = [pZ.G(:,temp == 0), 0.5*Gquad, pZ.Grest];

    % generate zonotope
    Z = zonotope([c,G]);
    
else
    
    Z = zonotope([pZ.c,pZ.Grest]);
    
end

%------------- END OF CODE --------------