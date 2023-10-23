function fs = generateRandom(varargin)
% generateRandom - generates a random full-dimensional space
%
% Syntax:
%    fs = generateRandom()
%    fs = generateRandom('Dimension',n)
%
% Inputs:
%    Name-Value pairs (all options, arbitrary order):
%       <'Dimension',n> - dimension
%
% Outputs:
%    fs - random fullspace
%
% Example: 
%    fs = fullspace.generateRandom();
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       25-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% name-value pairs -> number of input arguments is always a multiple of 2
if mod(nargin,2) ~= 0
    throw(CORAerror('CORA:evenNumberInputArgs'));
else
    % read input arguments
    NVpairs = varargin(1:end);
    % check list of name-value pairs
    checkNameValuePairs(NVpairs,{'Dimension'});
    % dimension given?
    [NVpairs,n] = readNameValuePair(NVpairs,'Dimension');
end

% default computation for dimension
if isempty(n)
    nmax = 10;
    n = randi(nmax);
end

% init fullspace object
fs = fullspace(n);

% ------------------------------ END OF CODE ------------------------------
