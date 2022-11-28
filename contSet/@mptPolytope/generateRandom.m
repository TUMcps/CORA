function P = generateRandom(varargin)
% generateRandom - Generates a random zonotope
%
% Syntax:  
%    poly = mptPolytope.generateRandom()
%    poly = mptPolytope.generateRandom('Dimension',n)
%    poly = mptPolytope.generateRandom('Dimension',n,'Center',c)
%
% Inputs:
%    Name-Value pairs (all options, arbitrary order):
%       <'Dimension',n> - dimension
%       <'Center',c> - center
%
% Outputs:
%    P - random polytope
%
% Example: 
%    P = mptPolytope.generateRandom('Dimension',2);
%    plot(P);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/generateRandom

% Author:       Niklas Kochdumper
% Written:      05-May-2020
% Last update:  19-May-2022 (MW, name-value pair syntax)
% Last revision:---

%------------- BEGIN CODE --------------

% name-value pairs -> number of input arguments is always a multiple of 2
if mod(nargin,2) ~= 0
    throw(CORAerror('CORA:evenNumberInputArgs'));
else
    % read input arguments
    NVpairs = varargin(1:end);
    % check list of name-value pairs
    checkNameValuePairs(NVpairs,{'Dimension','Center'});
    % dimension given?
    [NVpairs,n] = readNameValuePair(NVpairs,'Dimension');
    % center given?
    [NVpairs,c] = readNameValuePair(NVpairs,'Center');
end

% default dimension
if isempty(n)
    if isempty(c)
        maxdim = 10;
        n = randi(maxdim);
    else
        n = length(c);
    end
end

% default computation of center (only approximately center of P)
if isempty(c)
    c = -10 + 20*rand(n,1);
end

% generate random normal distribution
A = rand(n,n);
% compute random points
sigma = 0.5*(A+A') + n*eye(n);
nrPoints = n*10;
points = mvnrnd(c,sigma,nrPoints)';


% instantiate polytope
P = mptPolytope.enclosePoints(points);

%------------- END OF CODE --------------