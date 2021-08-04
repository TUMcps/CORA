function pZ = polyZonotope(obj)
% polyZonotope - convert a zonotope object to a polyZonotope object
%
% Syntax:  
%    pZ = polyZonotope(obj)
%
% Inputs:
%    obj - zonotope object (class zonotope)
%
% Outputs:
%    pZ - polynomial zonotope object (class polyZonotope)
%
% Example: 
%    zono = zonotope([1 2 0 -1;3 1 2 2]);
%    pZ = polyZonotope(zono);
%    
%    figure;
%    plot(zono,[1,2],'b','Filled',true,'EdgeColor','none');
%    xlim([-3,5]); ylim([-3,9]);
%
%    figure;
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%    xlim([-3,5]); ylim([-3,9]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval/polyZonotope, taylm/polyZonotope

% Author:       Niklas Kochdumper
% Written:      25-June-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    c = obj.Z(:,1);
    G = obj.Z(:,2:end);
    expMat = eye(size(G,2));

    pZ = polyZonotope(c,G,[],expMat);
    
%------------- END OF CODE --------------