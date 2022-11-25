function cZ = generateRandom(varargin)
% generateRandom - Generates a random constrained zonotope
%
% Syntax:  
%    cZ = generateRandom()
%    cZ = generateRandom(dim)
%    cZ = generateRandom(dim,cen)
%    cZ = generateRandom(dim,cen,nrOfGens)
%    cZ = generateRandom(dim,cen,nrOfGens,nrOfCons)
%
% Inputs:
%    dim      - dimension
%    cen      - zonotope center
%    nrOfGens - number of generators
%    nrOfCons - number of constraints
%
% Outputs:
%    cZ - random conZonotope object
%
% Example: 
%    cZ = conZonotope.generateRandom(2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/generateRandom

% Author:       Niklas Kochdumper
% Written:      30-Oct-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    dim = []; cen = []; nrOfGens = []; nrOfCons = [];

    if nargin >= 1 && ~isempty(varargin{1})
        dim = varargin{1};
    end
    if nargin >= 2 && ~isempty(varargin{2})
        cen = varargin{2};
    end
    if nargin >= 3 && ~isempty(varargin{3})
        nrOfGens = varargin{3};
    end
    if nargin >= 4 && ~isempty(varargin{4})
        nrOfCons = varargin{4};
    end
    
    % generate random zonotope
    Z = zonotope.generateRandom(dim,cen,nrOfGens);
    
    % get dimension and number of generators
    G = generators(Z);
    [n,m] = size(G);
    
    % select number of constraints
    if isempty(nrOfCons)
              
        nrOfCons = 0;
        
        if m > n
            nrOfCons = randi([0,m-n]);
        end
    else
        nrOfCons = min(nrOfCons,m-n);
    end
    
    % generate random constraints
    if nrOfCons > 0
        
        p = -1 + 2*rand(m,1);
        Aeq = -5 + 10*rand(nrOfCons,m);
        beq = Aeq*p;
        
        % construct constrained zonotope
        cZ = conZonotope(center(Z),G,Aeq,beq);       
        
    else
        cZ = conZonotope(Z); 
    end
end

%------------- END OF CODE --------------