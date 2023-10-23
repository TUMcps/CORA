function Z = generateRandom(varargin)
% generateRandom - Generates a random zonotope
%
% Syntax:
%    Z = zonotope.generateRandom()
%    Z = zonotope.generateRandom('Dimension',n)
%    Z = zonotope.generateRandom('Dimension',n,'NrGenerators',nrGens)
%
% Inputs:
%    Name-Value pairs (all options, arbitrary order):
%       <'Dimension',n> - dimension
%       <'Center',c> - center
%       <'NrGenerators',nrGens> - number of generators
%       <'Distribution',type> - distribution for generators
%           typeDist has to be {'uniform', 'exp', 'gamma'}
%
% Outputs:
%    Z - random zonotope
%
% Example: 
%    Z1 = zonotope.generateRandom();
%    Z2 = zonotope.generateRandom('Dimension',3);
%    Z3 = zonotope.generateRandom('Center',ones(2,1));
%    Z4 = zonotope.generateRandom('Dimension',4,'NrGenerators',10);
%    Z5 = zonotope.generateRandom('Distribution','gamma');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       17-September-2019
% Last update:   24-September-2019
%                01-May-2020 (integration of randomZonotope.m)
%                19-May-2022 (name-value pair syntax)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% name-value pairs -> number of input arguments is always a multiple of 2
if mod(nargin,2) ~= 0
    throw(CORAerror('CORA:evenNumberInputArgs'));
else
    % read input arguments
    NVpairs = varargin(1:end);
    % check list of name-value pairs
    checkNameValuePairs(NVpairs,{'Dimension','Center','NrGenerators','Distribution'});
    % dimension given?
    [NVpairs,n] = readNameValuePair(NVpairs,'Dimension');
    % center given?
    [NVpairs,c] = readNameValuePair(NVpairs,'Center');
    % number of generators given?
    [NVpairs,nrGens] = readNameValuePair(NVpairs,'NrGenerators');
    % distribution for generators given?
    [NVpairs,type] = readNameValuePair(NVpairs,'Distribution');
end


% default computation for dimension
if isempty(n)
    if isempty(c)
        nmax = 10;
        n = randi(nmax);
    else
        n = length(c);
    end
end

% default computation for center
if isempty(c)
    c = 10*randn(n,1);
end

% default number of generators
if isempty(nrGens)
    nrGens = 2*n;
end

% default distribution
if isempty(type)
    type = 'uniform';
end

% generate random vector for the generator lengths
l = zeros(nrGens,1);

% uniform distribution
if strcmp(type,'uniform')
    l = rand(nrGens,1);

% exponential distribution    
elseif strcmp(type,'exp')
    l = exprnd(1,nrGens,1);

% gaussian distribution    
elseif strcmp(type,'gamma')
    l = gamrnd(2,1,nrGens,1);
end

% init generator matrix
G = zeros(n,nrGens);
% create generators
for i=1:nrGens
    % generate random point on sphere
    gTmp = randomPointOnSphere(n);
    % stretch by length
    G(:,i) = l(i)*gTmp;
end

% instantiate zonotope
Z = zonotope(c,G);

end

% ------------------------------ END OF CODE ------------------------------
