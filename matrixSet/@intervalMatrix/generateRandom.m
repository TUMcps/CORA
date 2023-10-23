function intMat = generateRandom(varargin)
% generateRandom - generates a random interval matrix
%
% Syntax:
%    intMat = intervalMatrix.generateRandom()
%    intMat = intervalMatrix.generateRandom('Dimension',[n,m])
%    intMat = intervalMatrix.generateRandom('Dimension',[n,m],'Center',c)
%    intMat = intervalMatrix.generateRandom('Dimension',n,'Center',c,'MaxRadius',r)
%
% Inputs:
%    Name-Value pairs (all options, arbitrary order):
%       <'Dimension',[n,m]> - dimension
%       <'Center',c> - center
%       <'MaxRadius',r> - maximum radius for each dimension
%
% Outputs:
%    intMat - random zonotope
%
% Example: 
%    intMat1 = intervalMatrix.generateRandom();
%    intMat2 = intervalMatrix.generateRandom('Dimension',3);
%    intMat3 = intervalMatrix.generateRandom('Center',ones(2,1));
%    intMat4 = intervalMatrix.generateRandom('MaxRadius',ones(2,1));
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       03-April-2023
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
    checkNameValuePairs(NVpairs,{'Dimension','Center','MaxRadius'});
    % dimension given?
    [NVpairs,n] = readNameValuePair(NVpairs,'Dimension');
    % center given?
    [NVpairs,c] = readNameValuePair(NVpairs,'Center');
    % maximum radius given?
    [NVpairs,r] = readNameValuePair(NVpairs,'MaxRadius',@(val)all(all(val>=0)));
end


% default computation for dimension
if isempty(n)
    if isempty(c) && isempty(r)
        nmax = 10;
        n = randi(nmax,1,2);
    elseif ~isempty(c)
        n = size(c);
    elseif ~isempty(r)
        n = size(r);
    end
end

% default computation for center
if isempty(c)
    c = 10*randn(n);
end

% default computation of maximum radius
if isempty(r)
    r = 10*rand(n);
elseif isscalar(r)
    r = repmat(r,n(1),n(2));
end

% instantiate interval matrix
intMat = intervalMatrix(c,r);

% ------------------------------ END OF CODE ------------------------------
