function hyp = generateRandom(varargin)
% generateRandom - Generates a random conHyperplane
%
% Syntax:
%    hs = conHyperplane.generateRandom()
%    hs = conHyperplane.generateRandom('Dimension',n)
%    hs = conHyperplane.generateRandom('Dimension',n,'NormalVector',c)
%    hs = conHyperplane.generateRandom('Dimension',n,'NormalVector',c,'Offset',d)
%
% Inputs:
%    Name-Value pairs (all options, arbitrary order):
%       <'Dimension',n> - dimension
%       <'NormalVector',c> - normal vector of conHyperplane
%       <'Offset',d> - offset of conHyperplane
%       <'ContraintMatrix',C> - constraint matrix
%       <'ConstraintVector',d> - constraint vector
%
% Outputs:
%    hs - random conHyperplane
%
% Example: 
%    hyp1 = conHyperplane.generateRandom();
%    hyp2 = conHyperplane.generateRandom('Dimension',3);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: halfspace/generateRandom

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
    checkNameValuePairs(NVpairs,{'Dimension','NormalVector','Offset'});
    % dimension given?
    [NVpairs,n] = readNameValuePair(NVpairs,'Dimension');
    % normal vector given?
    [NVpairs,a] = readNameValuePair(NVpairs,'NormalVector');
    % offset given?
    [NVpairs,b] = readNameValuePair(NVpairs,'Offset');
    % constraint matrix given?
    [NVpairs,C] = readNameValuePair(NVpairs,'ContraintMatrix');
    % constraint vector given?
    [NVpairs,d] = readNameValuePair(NVpairs,'ConstraintVector');
end

% check input arguments
inputArgsCheck({{n,'att','numeric','nonnan'};
                {a,'att','numeric','nonnan'};
                {b,'att','numeric','nonnan'};
                {C,'att','numeric','nonnan'};
                {d,'att','numeric','nonnan'}});

% default computation for dimension
hs = halfspace.generateRandom('Dimension',n,'NormalVector',a,'Offset',b);
n = dim(hs);

% default computation for constraint matrix
if isempty(C)
    if isempty(d)
        m = n;
    else
        m = length(d);
    end

    C = 5 * randn(m,n);
    % normalize
    C = C ./ vecnorm(C,2);
end

% default computation for constraint vector
if isempty(d)
    d = 5*randn(size(C,1),1);
end

% instantiate interval
hyp = conHyperplane(hs.c',hs.d,C,d);

% ------------------------------ END OF CODE ------------------------------
