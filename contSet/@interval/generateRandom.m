function I = generateRandom(varargin)
% generateRandom - Generates a random interval
%
% Syntax:
%    I = interval.generateRandom()
%    I = interval.generateRandom('Dimension',n)
%    I = interval.generateRandom('Dimension',n,'Center',c)
%    I = interval.generateRandom('Dimension',n,'Center',c,'MaxRadius',r)
%
% Inputs:
%    Name-Value pairs (all options, arbitrary order):
%       <'Dimension',n> - dimension
%       <'Center',c> - center
%       <'MaxRadius',r> - maximum radius for each dimension or scalar
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

% Authors:       Mark Wetzlinger, Tobias Ladner
% Written:       17-September-2019
% Last update:   19-May-2022 (MW, name-value pair syntax)
%                23-February-2023 (MW, add 'Center' and 'MaxRadius')
%                22-May-2023 (TL, bugfix: all dimensions had same radius)
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
    [NVpairs,r] = readNameValuePair(NVpairs,'MaxRadius');
end

% check if the center matches the dimension (if both provided)
if ~isempty(n)
    if isscalar(n)
        % rewrite as [n,1] for easier handling of matrices
        n = [n,1];
    end
end

% default computation of dimension
if isempty(n)
    if ~isempty(c)
        % center given -> read out dimension
        n = size(c);
    elseif ~isempty(r) && ~isscalar(r)
        % radius given with dimension specs
        n = size(r);
    else
        nmax = 10;
        n = [randi(nmax),1];
    end
end

if CHECKS_ENABLED
    if ~isempty(c) && any(n ~= size(c))
        throw(CORAerror('CORA:wrongValue','name-value pair Center',...
            'has to match the dimension of name-value pair Dimension'));
    end
    if ~isempty(r) && ~isscalar(r) && any(n ~= size(r))
        throw(CORAerror('CORA:wrongValue','name-value pair MaxRadius',...
            'has to match the dimension of name-value pair Dimension'));
    end
end

% default computation of center
if isempty(c)
    % set somewhere in the neighborhood of the origin
    c = -2 + 4*rand(n);
end

% default computation of maximum radius
if isempty(r)
    r = 10*rand(n);
end

% effective radius of interval (has to be symmetric to maintain center)
rad = r/2 .* rand(n);

% instantiate interval
I = interval(c-rad,c+rad);

% ------------------------------ END OF CODE ------------------------------
