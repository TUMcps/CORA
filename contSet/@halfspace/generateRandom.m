function hs = generateRandom(varargin)
% generateRandom - Generates a random halfspace
%
% Syntax:
%    hs = halfspace.generateRandom()
%    hs = halfspace.generateRandom('Dimension',n)
%    hs = halfspace.generateRandom('Dimension',n,'NormalVector',c)
%    hs = halfspace.generateRandom('Dimension',n,'NormalVector',c,'Offset',d)
%
% Inputs:
%    Name-Value pairs (all options, arbitrary order):
%       <'Dimension',n> - dimension
%       <'NormalVector',c> - normal vector of halfspace
%       <'Offset',d> - offset of halfspace
%
% Outputs:
%    hs - random halfspace
%
% Example: 
%    hs1 = halfspace.generateRandom();
%    hs2 = halfspace.generateRandom('Dimension',3);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

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
    checkNameValuePairs(NVpairs,{'Dimension','NormalVector','Offset'});
    % dimension given?
    [NVpairs,n] = readNameValuePair(NVpairs,'Dimension');
    % normal vector given?
    [NVpairs,c] = readNameValuePair(NVpairs,'NormalVector');
    % offset given?
    [NVpairs,d] = readNameValuePair(NVpairs,'Offset');
end

% check input arguments
inputArgsCheck({{n,'att','numeric','nonnan'};
                {c,'att','numeric','nonnan'};
                {d,'att','numeric','nonnan'}});

% default computation for dimension
if isempty(n)
    if isempty(c)
        nmax = 10;
        n = randi(nmax);
    else
        n = length(c);
    end
end

% default computation for normal vector
if isempty(c)
    c = 5 * randn(n,1);
    % normalize
    c = c ./ vecnorm(c,2);
end

% default computation for offset
if isempty(d)
    d = 5*randn(1);
end


% instantiate interval
hs = halfspace(c,d);

% ------------------------------ END OF CODE ------------------------------
