function [E] = generateRandom(varargin)
% generateRandom - Generates a random ellipsoid
%
% Syntax:  
%    E = generateRandom() Generates an ellipsoid with random orientation,
%        center and dimension
%    E = generateRandom(isdegenerate) Generates a random
%        degenerate/non-degenerate ellipsoid based on given boolean
%        isdegenerate
%    E = generateRandom(isdegenerate,dim) Generates a random
%        degenerate/non-degenerate ellipsoid with dimension dim
%    E = generateRandom(dim,isdegenerate) see above
%
% Inputs:
%   (opt.) isdegenerate - type of ellipsoid
%
% Outputs:
%    E - random ellipsoid
%
% Example: 
%    E = ellipsoid.generateRandom(true);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann, Matthias Althoff
% Written:      13-March-2019
% Last update:  02-Sep-2019 (rename generate -> generateRandom)
%               19-March-2021 (complete rewrite)
%               30-July-2021 (removed "makedist" (to remove toolbox dep.))
% Last revision:---
%------------- BEGIN CODE --------------
isdegenerate=[];
maxdim = 30;
% choose random state dimension (if not specified by function call)
n = randi(maxdim);
% check which calling syntax is used
if nargin==1
    if isa(varargin{1},'logical')
        isdegenerate = varargin{1};
    elseif floor(varargin{1})==varargin{1} && varargin{1}>0
        n = varargin{1};
    end
elseif nargin==2
    if isa(varargin{1},'logical') && floor(varargin{2})==varargin{2} && varargin{2}>0
        isdegenerate = varargin{1};
        n = varargin{2};
    elseif isa(varargin{2},'logical') && floor(varargin{1})==varargin{1} && varargin{1}>0
        isdegenerate = varargin{2};
        n = varargin{1};
    else
        error('Wrong input arguments');
    end
elseif nargin>2
    error('Too many input arguments');
end
% generate a n-by-n random matrix (normal distribution) 
tmp = randn(n);
Q = tmp'*tmp;
% make sure Q is positive-semidefinite (i.e. a valid ellipsoid shape
% matrix)
Q = 1/2*(Q+Q');
E = ellipsoid;
TOL = E.TOL;
% if user explicitely specified degeneracy, keep going; otherwise good
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

% choose random center 
q = randn(n,1);
% construct resulting ellipsoid
E = ellipsoid(Q,q);
%------------- END OF CODE --------------