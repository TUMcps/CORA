function C = generateRandom(varargin)
% generateRandom - Generates a random capsule
%
% Syntax:
%    C = capsule.generateRandom()
%    C = capsule.generateRandom('Dimension',n)
%    C = capsule.generateRandom('Dimension',n,'Center',c)
%    C = capsule.generateRandom('Dimension',n,'Center',c,'Radius',r)
%
% Inputs:
%    Name-Value pairs (all options, arbitrary order):
%       <'Dimension',n> - dimension
%       <'Center',c> - center
%       <'Radius',r> - radius
%
% Outputs:
%    C - random capsule
%
% Example: 
%    C = capsule.generateRandom('Dimension',3);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       17-September-2019
% Last update:   19-May-2022 (name-value pairs syntax)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% name-value pairs -> number of input arguments is always a multiple of 2
if mod(nargin,2) ~= 0
    throw(CORAerror('CORA:evenNumberInputArgs'));
else
    % read input arguments
    NVpairs = varargin(1:end);
    % check list of name-value pairs
    checkNameValuePairs(NVpairs,{'Dimension','Center','Radius'});
    % dimension given?
    [NVpairs,n] = readNameValuePair(NVpairs,'Dimension');
    % center given?
    [NVpairs,c] = readNameValuePair(NVpairs,'Center');
    % radius given?
    [NVpairs,r] = readNameValuePair(NVpairs,'Radius');
end

% check input arguments
inputArgsCheck({{n,'att','numeric','nonnan'};
                {c,'att','numeric','nonnan'};
                {r,'att','numeric','nonnan'}});

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

% default computation for radius
if isempty(r)
    r_max = 5;
    r = rand * r_max;
end

% default computation for generator
g = rand(n,1);

% instantiate capsule
C = capsule(c,g,r);


% ------------------------------ END OF CODE ------------------------------
