function cPZ = quadMap(varargin)
% quadMap - computes the quadratic map of a constrained polynomial zonotope
%
% Syntax:
%    cPZ = quadMap(cPZ1,Q)
%    cPZ = quadMap(cPZ1,cPZ2,Q)
%
% Inputs:
%    cPZ1,cPZ2 - conPolyZono objects
%    Q - quadratic coefficients as a cell of matrices
%
% Outputs:
%    cPZ - conPolyZono object
%
% Example: 
%    c = [0;0];
%    G = [2 1 0;2 0 1];
%    E = [1 1 1; 0 1 0; 0 0 1; 0 0 0];
%    A = [1 1 -0.5];
%    b = 0.5;
%    EC = [0 0 0; 2 0 0; 0 2 0; 0 0 1];
%    cPZ = conPolyZono(c,G,E,A,b,EC);
%    
%    Q{1} = [1 2; -1 2];
%    Q{2} = [-3 0; 1 1];
%
%    res = quadMap(cPZ,Q);
%
%    figure; hold on
%    plot(cPZ,[1,2],'FaceColor','b','Splits',12);
%
%    figure; hold on
%    plot(res,[1,2],'FaceColor','r','Splits',12);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/quadMap, zonotope/quadMap, cubMap

% Authors:       Niklas Kochdumper
% Written:       19-January-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin == 1
    throw(CORAerror('CORA:notEnoughInputArgs',2));

elseif nargin == 2                                  % quadratic map
    % syntax:
    % cPZ = quadMap(cPZ1,Q)
    
    % compute quadratic map for polynomial zonotopes
    cPZ = varargin{1};
    pZ = polyZonotope(cPZ.c,cPZ.G,cPZ.GI,cPZ.E,cPZ.id);
    
    pZ = quadMap(pZ,varargin{2});
    
    % convert to constraint polynomial zonotope
    cPZ = conPolyZono(pZ.c,pZ.G,pZ.E,cPZ.A,cPZ.b,cPZ.EC,pZ.GI,pZ.id);
                  
elseif nargin == 3                                  % mixed quadratic map
    % syntax:
    % cPZ = quadMap(cPZ1,cPZ2,Q)
    
    % compute quadratic map for polynomial zonotopes
    cPZ1 = varargin{1}; cPZ2 = varargin{2};
    pZ1 = polyZonotope(cPZ1.c,cPZ1.G,cPZ1.GI,cPZ1.E,cPZ1.id);
    pZ2 = polyZonotope(cPZ2.c,cPZ2.G,cPZ2.GI,cPZ2.E,cPZ2.id);
    
    pZ = quadMap(pZ1,pZ2,varargin{3});
    
    % convert to constraint polynomial zonotope
    cPZ = conPolyZono(pZ);
    
    % update constraints
    cPZ = updateConstraints(cPZ,cPZ1,cPZ2);
    
else
    
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% ------------------------------ END OF CODE ------------------------------
