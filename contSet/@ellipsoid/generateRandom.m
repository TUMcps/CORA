function E = generateRandom(varargin)
% generateRandom - Generates a random ellipsoid
%
% Syntax:
%    E = ellipsoid.generateRandom()
%    E = ellipsoid.generateRandom('Dimension',n)
%    E = ellipsoid.generateRandom('Center',c)
%    E = ellipsoid.generateRandom('IsDegenerate',true)
%
% Inputs:
%    Name-Value pairs (all options, arbitrary order):
%       <'Dimension',n> - dimension
%       <'Center',c> - center
%       <'IsDegenerate',isdegenerate> - degeneracy
%
% Outputs:
%    E - random ellipsoid
%
% Example: 
%    E = ellipsoid.generateRandom('Dimension',2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann, Matthias Althoff
% Written:       13-March-2019
% Last update:   02-September-2019 (rename generate -> generateRandom)
%                19-March-2021 (complete rewrite)
%                30-July-2021 (removed "makedist" (to remove toolbox dep.))
%                19-May-2022 (MW, name-value pair syntax)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% name-value pairs -> number of input arguments is always a multiple of 2
if mod(nargin,2) ~= 0
    throw(CORAerror('CORA:evenNumberInputArgs'));
else
    % read input arguments
    NVpairs = varargin(1:end);
    % check list of name-value pairs
    checkNameValuePairs(NVpairs,{'Dimension','IsDegenerate','Center'});
    % dimension given?
    [NVpairs,n] = readNameValuePair(NVpairs,'Dimension');
    % degeneracy determined?
    [NVpairs,isdegenerate] = readNameValuePair(NVpairs,'IsDegenerate');
    % center given?
    [NVpairs,q] = readNameValuePair(NVpairs,'Center');
end

% default computation for dimension
if isempty(n)
    if isempty(q)
        maxdim = 30;
        n = randi(maxdim);
    else
        n = length(q);
    end
end

% default computation for center
if isempty(q)
    q = randn(n,1);
end

% default value for isdegenerate
if isempty(isdegenerate)
    isdegenerate = false;
end

% generate a n-by-n random matrix (normal distribution) 
tmp = randn(n);
Q = tmp'*tmp;
% make sure Q is positive-semidefinite -> valid ellipsoid shape matrix
Q = 1/2*(Q+Q');
E = ellipsoid.empty(1);
TOL = E.TOL;
% if user explicitly specified degeneracy, keep going; otherwise good
% enough
if ~isempty(isdegenerate)
    % to get a truly non-degenerate ellipsoid, we make sure to not
    % choose any singular values near 1/cond(Q)=TOL
    [U,S,~] = svd(Q);
    s = diag(S);
    % minimum singular value we allow for non-degenerate ellipsoids
    sv_min = sqrt(TOL)*max(s);
    max_sv = max(s);
    % too small for comfort
    ind_ts = s<sv_min;
    s(ind_ts) = max_sv + 1/2*(max_sv-sv_min)*rand(sum(ind_ts),1);
    
    % if should be degenerate
    if isdegenerate
        % choose random number of degenerate dimension (keep at least one
        % non-degenerate)
        n_d = randi(n);
        dims = randi(n_d,n_d,1);
        s(dims) = 0;
    end
    Q = U*diag(s)*U';
end

% make sure shape matrix is symmetric
Q = 1/2*(Q+Q');

% instantiate ellipsoid
E = ellipsoid(Q,q);

% ------------------------------ END OF CODE ------------------------------
