function I = generateRandom(varargin)
% generateRandom - Generates a random interval
%
% Syntax:  
%    I = interval.generateRandom()
%    I = interval.generateRandom('Dimension',n)
%
% Inputs:
%    Name-Value pairs (all options, arbitrary order):
%       <'Dimension',d> - dimension
%
% Outputs:
%    I - random interval
%
% Example: 
%    I = interval.generateRandom('Dimension',2);
%    plot(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      17-Sep-2019
% Last update:  19-May-2022 (name-value pair syntax)
% Last revision:---

%------------- BEGIN CODE --------------

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

% check input arguments
inputArgsCheck({{n,'att','numeric','nonnan'}});

% default computation for dimension
if isempty(n)
    nmax = 10;
    n = randi(nmax);
end

% default computation of bounds
lb = -10; mid = 0; ub = 10;
infi = lb + rand(n,1)*(mid - lb);
supr = mid + rand(n,1)*(ub - mid);

% instantiate interval
I = interval(infi,supr);

%------------- END OF CODE --------------