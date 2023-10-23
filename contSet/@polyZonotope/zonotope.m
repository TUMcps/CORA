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
%    Z = zonotope(pZ);
%
%    figure; hold on;
%    plot(pZ,[1,2],'Filled','r');
%    plot(Z,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, polytope

% Authors:       Niklas Kochdumper
% Written:       24-March-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if ~isempty(pZ.G)
    % determine dependent generators with exponents that are all even
    isEvenColumn = all(mod(pZ.E, 2) == 0, 1);
    Gquad = pZ.G(:, isEvenColumn);

    % compute zonotope parameter
    c = pZ.c + 0.5 * sum(Gquad,2);
    G = [pZ.G(:, ~isEvenColumn), 0.5*Gquad, pZ.GI];

    % generate zonotope
    Z = zonotope(c,G);
    
else
    % only independent generators
    Z = zonotope(pZ.c,pZ.GI);
    
end

% ------------------------------ END OF CODE ------------------------------
