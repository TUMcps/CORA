function zB = generateRandom(varargin)
% generateRandom - Generates a random zonotope bundle
%
% Syntax:
%    zB = zonoBundle.generateRandom()
%    zB = zonoBundle.generateRandom('Dimension',n)
%    zB = zonoBundle.generateRandom('Dimension',n,'NrZonotopes',nrZonos)
%
% Inputs:
%    Name-Value pairs (all options, arbitrary order):
%       <'Dimension',n> - dimension
%       <'NrZonotopes',nrZonos> - number of zonotopes in bundle
%
% Outputs:
%    zB - random zonotope bundle
%
% Example: 
%    zB = zonoBundle.generateRandom('Dimension',2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/generateRandom

% Authors:       Mark Wetzlinger, Tobias Ladner
% Written:       17-September-2019
% Last update:   19-May-2022 (MW, name-value pair syntax)
%                04-October-2024 (MW, faster algorithm)
%                22-May-2025 (TL, bug fix, fixed point was not always contained)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% name-value pairs -> number of input arguments is always a multiple of 2
if mod(nargin,2) ~= 0
    throw(CORAerror('CORA:evenNumberInputArgs'));
else
    % read input arguments
    NVpairs = varargin(1:end);
    % check list of name-value pairs
    checkNameValuePairs(NVpairs,{'Dimension','NrZonotopes'});
    % dimension given?
    [NVpairs,n] = readNameValuePair(NVpairs,'Dimension');
    % number of zonotopes given?
    [NVpairs,nrZonos] = readNameValuePair(NVpairs,'NrZonotopes');
end

% default computation for dimension
if isempty(n)
    nmax = 10;
    n = randi(nmax);
end

% default number of zonotopes
if isempty(nrZonos)
    nrZonos = 1+randi(4);
end

% algorithm for random generation of a non-empty zonotope bundle:
listZ = cell(nrZonos,1);
% 1) fix a point that will be contained in all generated zonotopes
p = randn(n,1);
for z=1:nrZonos
    % 2) randomly generate zonotopes centered at p
    Z = zonotope.generateRandom('Dimension',n,'Center',p);
    % 3) find vertex in random direction
    [~,v] = supportFunc(Z,rand(n,1));
    % 4) random factor [0,1] to shift the zonotope between v and p
    listZ{z} = Z + rand(1) * (v-p);
end

% instantiate zonotope bundle
zB = zonoBundle(listZ);

% ------------------------------ END OF CODE ------------------------------
