function res = isIntersecting_(cPZ,S,type,varargin)
% isIntersecting_ - determines if a constrained polynomial zonotope
%    intersects a set
%
% Syntax:
%    res = isIntersecting_(cPZ,S)
%    res = isIntersecting_(cPZ,S,type)
%
% Inputs:
%    cPZ - conPolyZono object
%    S - contSet object
%    type - type of check ('exact' or 'approx')
%
% Outputs:
%    res - true/false
%
% Example: 
%    c = [0;0];
%    G = [8 4; 0 8];
%    E = [1 0; 0 1; 0 0];
%    A = [1 1 -0.25];
%    b = 0.75;
%    EC = [2 0 0; 0 2 0; 0 0 1];
%    cPZ = conPolyZono(c,G,E,A,b,EC);
%    pZ1 = polyZonotope([0;0],0.5*[1 -2 1; 2 3 1],[],[1 0 2;0 1 1]);
%    pZ2 = polyZonotope([0;0],[1 -2 1; 2 3 1],[],[1 0 2;0 1 1]);
%
%    res1 = isIntersecting(cPZ,pZ1,'approx')
%    res2 = isIntersecting(cPZ,pZ2,'approx')
%
%    figure; hold on;
%    plot(cPZ,[1,2],'FaceColor','b','Splits',12);
%    plot(pZ1,[1,2],'g');
%    plot(pZ2,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/isIntersecting, polyZonotope/isIntersecting_, isempty

% Authors:       Niklas Kochdumper
% Written:       04-February-2021
% Last update:   ---
% Last revision: 27-March-2023 (MW, rename isIntersecting_)

% ------------------------------ BEGIN CODE -------------------------------

if strcmp(type,'exact')
    throw(CORAerror('CORA:noExactAlg',cPZ,S));
end

% call function for other set representations
if isa(S,'conPolyZono') || isa(S,'polyZonotope') || ...
   isa(S,'capsule')
    
    % fast check based on zonotope enclosure
    res = isIntersecting_(zonotope(cPZ),zonotope(S),type);
    
    if ~res
        return; 
    end
    
    % convert second set to constrained polynomial zonotope
    S = conPolyZono(S);
    
    % compute intersection of the two sets
    I = and_(cPZ,S,'exact');
    
    % check if the intersection is empty
    res = ~representsa_(I,'emptySet',eps);
    
    
elseif isa(S,'halfspace') || isa(S,'conHyperplane') || ...
       isa(S,'polytope') || isa(S,'zonotope') || ...
       isa(S,'interval') || isa(S,'zonoBundle') || ...
       isa(S,'conZonotope') || isa(S,'ellipsoid')

    res = isIntersecting_(S,cPZ,type);

else
    
    % throw error for given arguments
    throw(CORAerror('CORA:noops',cPZ,S));
end

% ------------------------------ END OF CODE ------------------------------
