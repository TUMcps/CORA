function cPZ = generateRandom(varargin)
% generateRandom - Generates a random constrained polynomial zonotope
%
% Syntax:
%    cPZ = conPolyZono.generateRandom()
%    cPZ = conPolyZono.generateRandom('Dimension',n)
%    cPZ = conPolyZono.generateRandom('Dimension',n,'nrGenerators',nrGens)
%    cPZ = conPolyZono.generateRandom('Dimension',n,'nrGenerators',nrGens,...
%           'NrFactors',nrFac)
%    cPZ = conPolyZono.generateRandom('Dimension',n,'nrGenerators',nrGens,...
%           'NrFactors',nrFac,'NrConstraints',nrCons)
%    cPZ = conPolyZono.generateRandom('Dimension',n,'nrGenerators',nrGens,...
%           'NrFactors',nrFac,'NrConstraints',nrCons,'NrIndGenerators',nrIndGens)
%
% Inputs:
%    Name-Value pairs (all options, arbitrary order):
%       <'Dimension',n> - dimension
%       <'NrGenerators',nrGens> - number of generators
%       <'NrConstraints',nrCons> - number of constraints
%       <'NrFactors',nrFac> - number of dependent factors
%       <'NrIndGenerators',nrIndGens> - number of independent generators
%
% Outputs:
%    cPZ - random constrained polynomial zonotope
%
% Example: 
%    cPZ = conPolyZono.generateRandom('Dimension',2,'NrGenerators',5)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/generateRandom

% Authors:       Niklas Kochdumper
% Written:       26-January-2020
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
    checkNameValuePairs(NVpairs,{'Dimension','NrGenerators',...
        'NrConstraints','NrFactors','NrIndGenerators'});
    % dimension given?
    [NVpairs,n] = readNameValuePair(NVpairs,'Dimension');
    % number of generators given?
    [NVpairs,nrGens] = readNameValuePair(NVpairs,'NrGenerators');
    % number of constraints given?
    [NVpairs,nrCons] = readNameValuePair(NVpairs,'NrConstraints');
    % number of factors given?
    [NVpairs,nrFac] = readNameValuePair(NVpairs,'NrFactors');
    % number of independent generators given?
    [NVpairs,nrIndGens] = readNameValuePair(NVpairs,'NrIndGenerators');
end

% check input arguments
inputArgsCheck({{n,'att','numeric','nonnan'};
                {nrGens,'att','numeric','nonnan'};
                {nrCons,'att','numeric','nonnan'};
                {nrFac,'att','numeric','nonnan'};
                {nrIndGens,'att','numeric','nonnan'}});
    
% generate random polynomial zonotope for the states
pZ = polyZonotope.generateRandom('Dimension',n,'NrGenerators',nrGens,...
    'NrFactors',nrFac,'NrIndGenerators',nrIndGens);

% determine number of constraints
nrFac = length(pZ.id);

if ~isempty(nrCons)
    nrCons = max(0,min(nrCons,nrFac-1)); 
else
    nrCons = randi([0 nrFac-1],1);
end

% generate random polynomial constraints 
if nrCons > 0
    pZcon = polyZonotope.generateRandom('NrGenerators',nrCons,...
        'NrFactors',nrFac);
else
    cPZ = conPolyZono(pZ); return;
end

% adapt the constraints to guarantee that the resulting set is non-empty
a = -1 + 2*rand(nrFac,1);
b = sum(pZcon.G.*prod(a.^pZcon.E,1),2);

% instantiate constrained polynomial zonotope
try
    cPZ = conPolyZono(pZ.c,pZ.G,pZ.E,pZcon.G,b,pZcon.E,pZ.GI);
catch
    cPZ = conPolyZono(pZ.c,pZ.G,pZ.E,pZ.GI);
end

% ------------------------------ END OF CODE ------------------------------
