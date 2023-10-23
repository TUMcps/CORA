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

% Authors:       Mark Wetzlinger
% Written:       17-September-2019
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

% default number of generators
if isempty(nrZonos)
    nrZonos = 1+randi(4);
end

% construct random zonotope bundle
listZ = cell(nrZonos,1);
listZ{1} = zonotope.generateRandom('Dimension',n);

% to ensure that zonotope bundle not empty, the center of the next zonotope
% in the list is a point contained in all previous zonotopes
for z=2:nrZonos
    inside = false;
    while ~inside
        inside = true;
        c = randPoint_(listZ{1},1,'standard');
        for i=2:z-1
            if ~contains_(listZ{i},c,'exact',0)
                inside = false; continue;
            end
        end
    end
    listZ{z,1} = zonotope.generateRandom('Dimension',n,'Center',c);
end

% instantiate zonotope bundle
zB = zonoBundle(listZ);

% ------------------------------ END OF CODE ------------------------------
