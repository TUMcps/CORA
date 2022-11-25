function res = restructureZono(pZ, order, method)
% restructureZono - Calculate a new representation of a polynomial zonotope
%                   through over-approximation with a linear zonotope
%
% Syntax:  
%    res = restructureZono(pZ, order)
%
% Inputs:
%    pZ - polyZonotope object
%    order - desired zonotope order of the dependent factors for the
%            resulting polynomial zonotope 
%    method - reduction technique for linear zonotopes
%             (see zonotope/reduce)
%
% Outputs:
%    res - resuling polyZonotope object which over-approximates pZ
%
% Example:
%    pZ = polyZonotope([0;0],[1 0 1;1 2 -2],[-1 0.1 -0.5;1.2 0.3 0.2],[1 0 1;0 1 2]);
%    pZnew1 = restructure(pZ,'zonotopeGirard',2);
%    pZnew2 = restructure(pZ,'zonotopeMethC',2);
%
%    hold on
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%    plot(zonotope(pZnew1),[1,2],'g');
%    plot(zonotope(pZnew2),[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/restructure

% Author:       Niklas Kochdumper
% Written:      25-July-2018 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % calculate zonotope over-approximation
    zono = zonotope(pZ);
    
    % reduce the zonotope to the desired order
    zono = reduce(zono,method,order);
    
    % construct the new polynomial zonotope object
    res = polyZonotope(zono);

%------------- END OF CODE --------------