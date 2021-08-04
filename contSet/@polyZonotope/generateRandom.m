function pZ = generateRandom(varargin)
% generateRandom - Generates a random polynomial zonotope
%
% Syntax:  
%    pZ = generateRandom()
%    pZ = generateRandom(dim)
%    pZ = generateRandom(dim,gen)
%    pZ = generateRandom(dim,gen,fac)
%    pZ = generateRandom(dim,gen,fac,ind)
%
% Inputs:
%    dim - (optional) dimension
%    gen - (optional) number of generators
%    fac - (optional) number of factors
%    ind - (optional) number of independent factors
%
% Outputs:
%    pZ - random polynomial zonotope
%
% Example: 
%    pZ = polyZonotope.generateRandom(2);
%    plot(pZ,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/generateRandom

% Author:       Niklas Kochdumper
% Written:      29-December-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % define bounds for random values
    dim_low = 2; dim_up = 10;
    orderGen_low = 1; orderGen_up = 10;
    fac_low = 1; fac_up = 5;
    gen_low = -2; gen_up = 2;
    exp_low = 1; exp_up = 4;

    % parse input arguments
    if nargin >= 1 && ~isempty(varargin{1})
        n = varargin{1};
    else
        n = dim_low + floor(rand(1) * (dim_up - dim_low + 1));
    end
    if nargin >= 2 && ~isempty(varargin{2})
        gen = varargin{2};
    else
        numGen_low = floor(orderGen_low*n);
        numGen_up = floor(orderGen_up*n);
        gen = numGen_low + floor(rand(1) * (numGen_up - numGen_low));
    end
    if nargin >= 3 && ~isempty(varargin{3})
        fac = varargin{3};
    else
        fac = floor(fac_low + rand(1)*(fac_up - fac_low));
    end
    if nargin >= 4 && ~isempty(varargin{4})
        indGens = varargin{4};
    else
        temp = randi([0,1]);
        if temp
            indGens = round(rand()*(orderGen_up*n));
        else
            indGens = 0;
        end
    end
    
    % fix number of generators
    if fac > gen
        gen = fac + 1;
    end

    % center vector
    c = gen_low + rand(n,1) * (gen_up - gen_low);

    % generator matrix
    G = gen_low + rand(n,gen) * (gen_up - gen_low);

    % exponent matrix
    expMat = zeros(fac,gen);
    
    for i = 1:round(fac*gen/2)
       
        % select random matrix entry
        ind1 = round(1 + rand()*(fac-1));
        ind2 = round(1 + rand()*(gen-1));
        
        % select random value
        expMat(ind1,ind2) = randomExponent(exp_low,exp_up);   
    end
    
    ind = find(sum(expMat,2) == 0);
    for i = 1:length(ind)
       temp = floor(1 + rand()*(size(expMat,2)-1));
       expMat(ind(i),temp) = randomExponent(exp_low,exp_up);
    end

    % remove redundant exponents
    [expMat,G] = removeRedundantExponents(expMat,G);
    
    ind = find(sum(expMat,1) == 0);
    
    if ~isempty(ind)
       expMat(:,ind) = [];
       G(:,ind) = [];
    end
    
    % generate random independent generators
    Grest = [];
    
    if indGens ~= 0
       Grest =  gen_low + rand(n,indGens) * (gen_up - gen_low);
    end

    % instantiate polyZonotope object
    pZ = polyZonotope(c,G,Grest,expMat);
    
end

% Auxiliary Functions -----------------------------------------------------

function e = randomExponent(l,u)
% generate a random exponent using a quadratic probability distribution
% since exponents have a large probability of being close to 1

    x = rand();
    e = l + (u-l)*(x-1)^2;
    e = round(e);
end

%------------- END OF CODE --------------