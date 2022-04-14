function cPZ = generateRandom(varargin)
% generateRandom - Generates a random constrained polynomial zonotope
%
% Syntax:  
%    cPZ = generateRandom()
%    cPZ = generateRandom(dim)
%    cPZ = generateRandom(dim,gen)
%    cPZ = generateRandom(dim,gen,fac)
%    cPZ = generateRandom(dim,gen,fac,cons)
%    cPZ = generateRandom(dim,gen,fac,cons,ind)
%
% Inputs:
%    dim - (optional) dimension
%    gen - (optional) number of generators
%    fac - (optional) number of factors
%    cons - (optional) number of constraints
%    ind - (optional) number of independent generators
%
% Outputs:
%    cPZ - random polynomial zonotope
%
% Example: 
%    cPZ = conPolyZono.generateRandom(2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/generateRandom

% Author:       Niklas Kochdumper
% Written:      26-January-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    n = []; gen = []; fac = []; cons = []; indGens = [];
    
    if nargin >= 1 && ~isempty(varargin{1})
        n = varargin{1};
    end
    if nargin >= 2 && ~isempty(varargin{2})
        gen = varargin{2};
    end
    if nargin >= 3 && ~isempty(varargin{3})
        fac = varargin{3};
    end
    if nargin >= 4 && ~isempty(varargin{4})
        cons = varargin{4}; 
    end
    if nargin >= 5 && ~isempty(varargin{5})
        indGens = varargin{5}; 
    end
    
    % generate random polynomial zonotope for the states
    pZ = polyZonotope.generateRandom(n,gen,fac,indGens);
    
    % determine number of constraints
    fac = length(pZ.id);
    
    if ~isempty(cons)
        cons = max(0,min(cons,fac-1)); 
    else
        cons = randi([0 fac-1],1);
    end
    
    % generate random polynomial constraints 
    if cons ~= 0
        pZcon = polyZonotope.generateRandom(cons,[],fac);
    else
        cPZ = conPolyZono(pZ); return;
    end
    
    % adapt the constraints such that it is guaranteed that the resulting
    % set is not empty
    a = -1 + 2*rand(fac,1);
    
    b = sum(pZcon.G.*prod(a.^pZcon.expMat,1),2);
    
    % construct resulting conPolyZono object
    try
        cPZ = conPolyZono(pZ.c,pZ.G,pZ.expMat,pZcon.G,b,pZcon.expMat,pZ.Grest);
    catch
        cPZ = conPolyZono(pZ.c,pZ.G,pZ.expMat,pZ.Grest);
    end   
end
    
%------------- END OF CODE --------------