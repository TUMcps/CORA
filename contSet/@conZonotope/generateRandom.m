function cZ = generateRandom(varargin)
% generateRandom - Generates a random constrained zonotope
%
% Syntax:
%    cZ = conZonotope.generateRandom()
%    cZ = conZonotope.generateRandom('Dimension',n)
%    cZ = conZonotope.generateRandom('Dimension',n,'NrGenerators',nrGens)
%    cZ = conZonotope.generateRandom('Dimension',n,'NrGenerators',nrGens,...
%           'NrConstraints',nrCons)
%
% Inputs:
%    Name-Value pairs (all options, arbitrary order):
%       <'Dimension',n> - dimension
%       <'Center',c> - center of the zonotope
%       <'NrGenerators',nrGens> - number of generators
%       <'NrConstraints',nrCons> - number of constraints
%
% Outputs:
%    cZ - random constrained zonotoope
%
% Example: 
%    cZ = conZonotope.generateRandom('Dimension',2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/generateRandom

% Authors:       Niklas Kochdumper
% Written:       30-October-2020
% Last update:   19-May-2022 (name-value pair syntax)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% name-value pairs -> number of input arguments is always a multiple of 2
if mod(nargin,2) ~= 0
    throw(CORAerror('CORA:evenNumberInputArgs'));
else
    % read input arguments
    NVpairs = varargin(1:end);
    % check list of name-value pairs
    checkNameValuePairs(NVpairs,{'Dimension','Center','NrGenerators','NrConstraints'});
    % dimension given?
    [NVpairs,n] = readNameValuePair(NVpairs,'Dimension');
    % center given?
    [NVpairs,c] = readNameValuePair(NVpairs,'Center');
    % number of generators given?
    [NVpairs,nrGens] = readNameValuePair(NVpairs,'NrGenerators');
    % number of constraints given?
    [NVpairs,nrCons] = readNameValuePair(NVpairs,'NrConstraints');
end

% check input arguments
inputArgsCheck({{n,'att','numeric','nonnan'};
                {c,'att','numeric','nonnan'};
                {nrGens,'att','numeric','nonnan'};
                {nrCons,'att','numeric','nonnan'}});

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
    c = 5*randn(n,1);
end

% default number of generators
if isempty(nrGens)
    nrGens = 2*n;
end

% default number of constraints
if isempty(nrCons)
    nrCons = 0;
    if nrGens > n
        nrCons = randi([0,nrGens-n]);
    end
end

% generate random zonotope
Z = zonotope.generateRandom('Dimension',n,'Center',c,'NrGenerators',nrGens);

% get generator matrix
G = generators(Z);

% generate random constraints
if nrCons > 0
    
    p = -1 + 2*rand(nrGens,1);
    Aeq = -5 + 10*rand(nrCons,nrGens);
    beq = Aeq*p;
    
    % construct constrained zonotope
    cZ = conZonotope(c,G,Aeq,beq);       
    
else
    cZ = conZonotope(Z); 
end

% ------------------------------ END OF CODE ------------------------------
