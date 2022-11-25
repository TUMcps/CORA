function pZ = polyZonotope(obj)
% polyZonotope - convert an interval object to a polynomial zonotope 
%
% Syntax:  
%    pZ = polyZonotope(obj)
%
% Inputs:
%    obj - interval object (class interval)
%
% Outputs:
%    pZ - polynomial zonotope object (class polyZonotope)
%
% Example: 
%    inter = interval([1;2],[3;5]);
%    pZ = polyZonotope(inter);
%    
%    plot(inter,[1,2],'b','Filled',true,'EdgeColor','none');
%    xlim([0,4]);
%    ylim([1,6]);
%
%    figure
%    plot(pZ, [1,2],'r','Filled',true,'EdgeColor','none');
%    xlim([0,4]);
%    ylim([1,6]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/polyZonotope, taylm/polyZonotope

% Author:       Niklas Kochdumper
% Written:      25-June-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

c = center(obj);
G = diag(rad(obj));
expMat = eye(length(c));

pZ = polyZonotope(c,G,[],expMat);
    
%------------- END OF CODE --------------