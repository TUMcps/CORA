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
%    cZ = conZonotope.generateRandom(2);
%    plot(cZ);
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
        
        options = optimoptions('linprog','display','off');
        
        % initialize valid factor domain
        l = -ones(m,1);
        u = ones(m,1);
        Aeq = []; beq = [];
        
        % loop over all constraints
        for i = 1:nrOfCons
           
            % select random direction
            d = rand(m,1) - 0.5*ones(m,1);
            
            % compute upper and lower bound of valid domain in direction
            [~,lb] = linprog(d,[],[],Aeq,beq,l,u,options);
            [~,ub] = linprog(-d,[],[],Aeq,beq,l,u,options);
            ub = -ub;
            
            % select random offset betwenn lower and upper bound
            b_ = lb + (ub-lb).*rand(1,1);
            
            % add constraint
            Aeq = [Aeq; d'];
            beq = [beq; b_];
        end
        
        % construct constrained zonotope
        cZ = conZonotope(center(Z),G,Aeq,beq);       
        
    else
        cZ = conZonotope(Z); 
    end
end

%------------- END OF CODE --------------