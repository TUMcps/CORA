function O = generateRandom(varargin)
% generateRandom - generates a random empty set
%
% Syntax:
%    O = generateRandom()
%    O = generateRandom('Dimension',n)
%
% Inputs:
%    Name-Value pairs (all options, arbitrary order):
%       <'Dimension',n> - dimension
%
% Outputs:
%    O - random emptySet
%
% Example: 
%    O = emptySet.generateRandom();
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner
% Written:       02-August-2023
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

% init emptySet object
O = emptySet(n);

% ------------------------------ END OF CODE ------------------------------
