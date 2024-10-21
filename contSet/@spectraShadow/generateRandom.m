function SpS = generateRandom(varargin)
% generateRandom - Generates a (somewhat biased) random spectrahedral
%    shadow.
%    Specifically, the spectrahedral shadow that is generated will be the
%    image under an affine map of the minkowski sum of an ellipsoid and a
%    polytope
%
% Syntax:
%    SpS = spectraShadow.generateRandom()
%    SpS = spectraShadow.generateRandom('Dimension',n)
%    SpS = spectraShadow.generateRandom('Dimension',n,'NrGenerators',nrGens)
%    SpS = spectraShadow.generateRandom('Dimension',n,'NrGenerators',nrGens,...
%           'IsBounded',isBounded)
%    SpS = spectraShadow.generateRandom('Dimension',n,'NrGenerators',nrGens,...
%           'IsBounded',isBounded,'IsDegenerate')
%
% Inputs:
%    Name-Value pairs (all options, arbitrary order):
%       <'Dimension',n> - dimension
%       <'Center',c> - center
%       <'NrGenerators',nrGens> - number of generators
%       <'IsBounded',isbounded> - boundedness (true/false)
%       <'IsDegenerate',isdegenerate> - degeneracy (true/false for
%           degenerate/non-degenerate. Non-degenerate is default)
%       <'Distribution',type> - distribution for generators
%           typeDist has to be {'uniform', 'exp', 'gamma'}
%
% Outputs:
%    SpS - random spectraShadow object
%
% Example: 
%    S1 = spectraShadow.generateRandom();
%    S2 = spectraShadow.generateRandom('Dimension',3);
%    S3 = spectraShadow.generateRandom('Center',ones(2,1));
%    S4 = spectraShadow.generateRandom('Dimension',4,'NrGenerators',10);
%    S5 = spectraShadow.generateRandom('Distribution','gamma');
%    S6 = spectraShadow.generateRandom('IsBounded',false);
%    S7 = spectraShadow.generateRandom('IsDegenerate',true);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Adrian Kulmburg
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
    checkNameValuePairs(NVpairs,{'Dimension','Center','NrGenerators','Distribution','IsBounded','IsDegenerate'});
    % dimension given?
    [NVpairs,n] = readNameValuePair(NVpairs,'Dimension',@(x) mod(x,1) == 0 && x > 0);
    % center given?
    [NVpairs,c] = readNameValuePair(NVpairs,'Center',@(v) isvector(v));
    % number of generators given?
    [NVpairs,nrGens] = readNameValuePair(NVpairs,'NrGenerators',@(x) mod(x,1) == 0 && x > 0);
    % distribution for generators given?
    [NVpairs,type] = readNameValuePair(NVpairs,'Distribution');
    % boundedness given?
    [NVpairs,isBnd] = readNameValuePair(NVpairs,'IsBounded','islogical');
    % degeneracy given?
    [NVpairs,isDeg] = readNameValuePair(NVpairs,'IsDegenerate','islogical');
end

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
    c = 10*randn(n,1);
end

% default number of generators
if isempty(nrGens)
    nrGens = 2*n;
end

% default distribution
if isempty(type)
    type = 'uniform';
end

% generate random vector for the generator lengths
l = zeros(nrGens,1);

% uniform distribution
if strcmp(type,'uniform')
    l = rand(nrGens,1);

% exponential distribution    
elseif strcmp(type,'exp')
    l = exprnd(1,nrGens,1);

% gaussian distribution    
elseif strcmp(type,'gamma')
    l = gamrnd(2,1,nrGens,1);
end

% default degeneracy
if isempty(isBnd)
    isBnd = true;
end

% default degeneracy
if isempty(isDeg)
    isDeg = false;
end

preimage_dim = floor(nrGens/2);

% init generator matrix
G = zeros([n preimage_dim]);
while rank(G) < min([n preimage_dim])
    % make sure G is full rank
    % create generators
    for i=1:preimage_dim
        % generate random point on sphere
        gTmp = randn([n 1]);
        gTmp = gTmp ./ norm(gTmp);
        % stretch by length
        G(:,i) = l(i)*gTmp;
    end
end

% Instantiate ellipsoid
E = ellipsoid.generateRandom('Dimension',preimage_dim);
% Instantiate polytope
P = polytope.generateRandom('Dimension',preimage_dim,'IsBounded',isBnd);

% If we need a degenerate spectrahedron
if isDeg && isBnd
    % Easy, we just select a random number of dimensions to delete
    r = rank(G);
    s = randi(r);

    [U, S, V] = svd(G);
    for i=1:s
        S((r-i+1),(r-i+1)) = 0;
    end
    G = U * S * V';
elseif isDeg && ~isBnd
    % Here it is harder, because we need to make sure that we keep the
    % unboundedness
    % So we will only 'degenerate' one direction
    [U, S, V] = svd(G);

    r = rank(G);

    for i=1:r
        S_copy = S;
        S_copy((r-i+1),(r-i+1)) = 0;

        G = U * S_copy * V';

        SpS_test = G * (spectraShadow(E) + spectraShadow(P)) + c;

        if ~isBounded(SpS_test)
            break % We found our degenerate, unbounded spectrahedron
        end

    end
    % If we didn't find it... well, tough luck, let's prioritize
    % unboundedness
    G = U * S * V';
end


% Combine everything we have to cook up some spectrahedral shadow
SpS = G * (spectraShadow(E) + spectraShadow(P)) + c;

% If nrGens is odd, add some random 1D zonotope
if ~even(nrGens)
    Z = zonotope(zeros([preimage_dim 1]), randn([preimage_dim 1]));
    SpS = SpS + G * Z;
end

end

% ------------------------------ END OF CODE ------------------------------
