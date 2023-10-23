function pZ = polyZonotope(I)
% polyZonotope - convert an interval object to a polynomial zonotope 
%
% Syntax:
%    pZ = polyZonotope(I)
%
% Inputs:
%    I - interval object (class interval)
%
% Outputs:
%    pZ - polynomial zonotope object (class polyZonotope)
%
% Example: 
%    I = interval([1;2],[3;5]);
%    pZ = polyZonotope(I);
%    
%    figure; hold on; xlim([0,4]); ylim([1,6]);
%    plot(I,[1,2],'FaceColor','b');
%    plot(pZ, [1,2],'FaceColor','r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/polyZonotope, taylm/polyZonotope

% Authors:       Niklas Kochdumper
% Written:       25-June-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

c = center(I);
G = diag(rad(I));
E = eye(length(c));

pZ = polyZonotope(c,G,[],E);
    
% ------------------------------ END OF CODE ------------------------------
