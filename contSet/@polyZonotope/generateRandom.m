function pZ = generateRandom(varargin)
% generateRandom - Generates a random polynomial zonotope
%
% Syntax:
%    pZ = polyZonotope.generateRandom()
%    pZ = polyZonotope.generateRandom('Dimension',n)
%    pZ = polyZonotope.generateRandom('Dimension',n,'NrGenerators',nrGens)
%    pZ = polyZonotope.generateRandom('Dimension',n,'NrGenerators',nrGens,...
%           'NrFactors',nrFac)
%    pZ = polyZonotope.generateRandom('Dimension',n,'NrGenerators',nrGens,...
%           'NrFactors',nrFac,'nrIndGenerators',nrIndGens)
%
% Inputs:
%    Name-Value pairs (all options, arbitrary order):
%       <'Dimension',n> - dimension
%       <'NrGenerators',nrGens> - number of generators
%       <'NrFactors',nrFac> - number of dependent factors
%       <'NrIndGenerators',nrIndGens> - number of independent generators
%
% Outputs:
%    pZ - random polynomial zonotope
%
% Example: 
%    pZ = polyZonotope.generateRandom('Dimension',2);
%    plot(pZ,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/generateRandom

% Authors:       Niklas Kochdumper
% Written:       29-December-2019
% Last update:   19-May-2022 (MW, name-value syntax)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% name-value pairs -> number of input arguments is always a multiple of 2
if mod(nargin,2) ~= 0
    throw(CORAerror('CORA:evenNumberInputArgs'));
else
    % read input arguments
    NVpairs = varargin(1:end);
    % check list of name-value pairs
    checkNameValuePairs(NVpairs,{'Dimension','NrGenerators','NrFactors','NrIndGenerators'});
    % dimension given?
    [NVpairs,n] = readNameValuePair(NVpairs,'Dimension');
    % number of generators given?
    [NVpairs,nrGens] = readNameValuePair(NVpairs,'NrGenerators');
    % number of factors given?
    [NVpairs,nrFac] = readNameValuePair(NVpairs,'NrFactors');
    % number of independent generators given?
    [NVpairs,nrIndGens] = readNameValuePair(NVpairs,'NrIndGenerators');
end

% check input arguments
inputArgsCheck({{n,'att','numeric','nonnan'};
                {nrGens,'att','numeric','nonnan'};
                {nrFac,'att','numeric','nonnan'};
                {nrIndGens,'att','numeric','nonnan'}});

% default computation for dimension
if isempty(n)
    nmax = 10;
    n = randi(nmax);
end

% default number of generators
if isempty(nrGens)
    nrGens = randi([3,10]);
end

% default number of factors
if isempty(nrFac)
    nrFac = Inf;
    while nrFac > nrGens + 1
        nrFac = randi([2,5]);
    end
end

% default number of independent generators
if isempty(nrGens)
    nrIndGens = randi([2,5]);
end

% define bounds for random values
gen_low = -2; gen_up = 2;
exp_low = 1; exp_up = 4;

% dependent generator matrix
G = gen_low + rand(n,nrGens) * (gen_up - gen_low);

% center vector
c = gen_low + rand(n,1) * (gen_up - gen_low);

% exponent matrix
E = zeros(nrFac,nrGens);

for j = 1:nrGens

    % select number of factors affecting j-th generator
    k = randi(nrFac);
    % select affected entries
    idxs = randperm(nrFac);
    idxs = idxs(1:k);

    for jj=1:k
        % select random value
        E(idxs(jj),j) = aux_randomExponent(exp_low,exp_up);   
    end
    
end

ind = find(sum(E,2) == 0);
for i = 1:length(ind)
   temp = floor(1 + rand*(size(E,2)-1));
   E(ind(i),temp) = aux_randomExponent(exp_low,exp_up);
end

% remove redundant exponents
[E,G] = removeRedundantExponents(E,G);

% generate random independent generators
GI = [];

if nrIndGens > 0
    GI = gen_low + rand(n,nrIndGens) * (gen_up - gen_low);
end

% instantiate polyZonotope object
pZ = polyZonotope(c,G,GI,E);
    

end


% Auxiliary functions -----------------------------------------------------

function e = aux_randomExponent(l,u)
% generate a random exponent using a quadratic probability distribution
% since exponents have a large probability of being close to 1

    x = rand();
    e = l + (u-l)*(x-1)^2;
    e = round(e);
end

% ------------------------------ END OF CODE ------------------------------
