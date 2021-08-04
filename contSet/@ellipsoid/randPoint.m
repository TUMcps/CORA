function p = randPoint(E,varargin)
% randPoint - generates a random point within an ellipsoid
%
% Syntax:  
%    p = randPoint(obj)
%    p = randPoint(obj,N)
%    p = randPoint(obj,type)
%    p = randPoint(obj,N,type)
%    p = randPoint(obj,N,'gaussian',pr)
%    p = randPoint(obj,'gaussian',pr)
%
% Inputs:
%    obj - ellipsoid object
%    N - number of random points
%    type - type of the random point ('standard' or 'gaussian' or 'extreme')
%    pr - probability that a value is within the set (only type = 'gaussian')
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    E = ellipsoid([1;0],rand(2,5));
%    p = randPoint(E);
% 
%    plot(E); hold on;
%    scatter(p(1,:),p(2,:),16,'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval/randPoint

% Author:        Victor Gassmann
% Written:       18-March-2021
% Last update:   25-June-2021 (MP, add type gaussian)
% Last revision: ---

%------------- BEGIN CODE --------------
% NOTE: This function does not produce a uniform distribution!
% parse input arguments
if ~isa(E,'ellipsoid')
    error('First argument has to be of type "ellipsoid"!');
end
defaultType = 'standard';
defaultN = 1;
types = {'standard','gaussian'};
if isempty(varargin)
    N = 1;
    type = defaultType;
elseif length(varargin)==1
    if isa(varargin{1},'double') && isscalar(varargin{1}) && mod(varargin{1},1)==0
        N = varargin{1};
        type = defaultType;
    elseif isa(varargin{1},'char') && any(strcmp(varargin{1},types))
        N = defaultN;
        type = varargin{1};
    else
        error('Wrong type for second input argument!');
    end
elseif length(varargin)==2
    if isa(varargin{1},'char') && strcmp(varargin{1},'gaussian') && isa(varargin{2},'double')
        type = varargin{1};
        pr = varargin{2};
    elseif ~(isa(varargin{1},'double') && isscalar(varargin{1}) && mod(varargin{1},1)==0)...
            || ~isa(varargin{2},'char') 
        error('Wrong type of input arguments!');
    else
        N = varargin{1};
        type = varargin{2};
    end
elseif length(varargin)==3
    if isnumeric(varargin{1}) && isscalar(varargin{1}) && isa(varargin{2},'char') ...
            && strcmp(varargin{2},'gaussian') && isa(varargin{3},'double')
        N = varargin{1};
        type = varargin{2};
        pr = varargin{3};
    else
        error('Type ''standard'' only supports up to 3 input arguments!');
    end
else
    error('Function only supports up to 4 input arguments!');
end


if isempty(E)
    p = [];
    return;
end
if rank(E)==0
    p = repmat(E.q,1,N);
    return;
end
c = E.q;
E = E + (-c);
r = rank(E);
n = dim(E);
n_rem = n-r;
[T,~,~] = svd(E.Q);
E = T'*E;
% if degenerate, project
E = project(E,1:r);
G = inv(sqrtm(E.Q));
E = G*E;



% generate different types of extreme points
if strcmp(type,'gaussian')
    if nargin == 3
        pt = randPoint@contSet(E,type,pr);
    else
        pt = randPoint@contSet(E,N,type,pr);
    end
    
    % stack again, backtransform and shift
    p = T*[inv(G)*pt;zeros(n_rem,N)] + c;
    
elseif strcmp(type,'standard')
    % generate points uniformely distributed (with N -> infinity)
    % on the unit hypersphere
    X = randn(dim(E),N);
    pt = X./sqrt(sum(X.^2,1));
    
    S = 2*rand(1,N)-1;
    pt = S.*pt;
    
    % stack again, backtransform and shift
    p = T*[inv(G)*pt;zeros(n_rem,N)] + c;
elseif strcmp(type,'extreme')
    pt = boundary(E,N);
    
    % stack again, backtransform and shift
    p = T*[inv(G)*pt;zeros(n_rem,N)] + c;
else
    error('"Type" has to be ''gaussian'',''standard'', or ''extreme''!');
end
%------------- END OF CODE --------------