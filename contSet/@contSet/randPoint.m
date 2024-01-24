function x = randPoint(S,varargin)
% randPoint - generates a random point within a given continuous set
%
% Syntax:
%    x = randPoint(S)
%    x = randPoint(S,N)
%    x = randPoint(S,N,'gaussian',pr)
%
% Inputs:
%    S - contSet object
%    N - number of random points
%    type - type of the random point ('standard', 'extreme', 'gaussian',
%           'uniform' or 'uniform:hitAndRun', 'uniform:billiardWalk')
%    pr - probability that a value is within the set (only type = 'gaussian') 
%
% Outputs:
%    x - random point in R^n
%
% Reference:
%    [1] Siotani, Minoru (1964). "Tolerance regions for a multivariate 
%        normal population". Annals of the Institute of Statistical 
%        Mathematics. 16 (1): 135–153. 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       19-November-2020
% Last update:   19-June-2021 (MP, generalization)
%                08-March-2023 (MW, reduce number of contains-calls)
%                28-March-2023 (MW, restructure relation to subclass)
%                22-May-2023 (AK, added 'uniform' sampling methods)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check number of input arguments
if nargin < 1
    throw(CORAerror('CORA:notEnoughInputArgs',1));
elseif nargin > 4
    throw(CORAerror('CORA:tooManyInputArgs',4));
end

% set default values for number of samples, method, and probability
[N,type,pr] = setDefaultValues({1,'standard',0.7},varargin);

% N can be numeric or 'all'
if isnumeric(N)
    checkN = {N,'att','numeric',{'scalar','integer','positive'}};
else
    checkN = {N,'str','all'};
end
% check input arguments
inputArgsCheck({{S,'att','contSet'};
                 checkN; % see above...
                {type,'str',{'standard','extreme','gaussian','uniform','uniform:hitAndRun','uniform:billiardWalk','radius'}};
                {pr,'att','numeric',{'<=',1,'>=',0}}});

% if N = 'all', then type has to be 'extreme' (nargin ensures that type has
% been set by calling function)
if ischar(N) && strcmp(N,'all') && nargin >= 3 && ~strcmp(type,'extreme')
    throw(CORAerror('CORA:wrongValue','third',...
        "If the number of points is 'all', the type has to be 'extreme'."));
end

% type = 'gaussian' implemented in contSet, other types in subclass methods
if ~strcmp(type,'gaussian')
    try
        x = randPoint_(S,N,type);
    catch ME
        if representsa_(S,'emptySet',eps)
            x = double.empty(dim(S),0);
        elseif representsa_(S,'origin',eps)
            x = repmat(zeros(dim(S),1),1,N);
        else
            rethrow(ME);
        end
    end
    return
end


% randPoint with type gaussian currently not supported for all contSet
% subclasses, therefore check if subclass is supported
if ~(isa(S,'zonotope') || isa(S,'interval') || isa(S,'ellipsoid') || isa(S,'polytope'))
    throw(CORAerror('CORA:notSupported',"The function randPoint for " + class(S) + ...
        " does not support the type = 'gaussian'."));
end

% zonotope: set does not have any generators
if isa(S,'zonotope') && isempty(generators(compact_(S,'zeros',eps)))
    x = S.c; return;
end

% generates a random vector according to Gaussian distribution within a
% given set enclose set by ellipsoid
if ~isa(S,'ellipsoid')
    E = ellipsoid(S);
else
    % degeneracy handling for ellipsoid: projection on proper subspace,
    % back-transformation after sampling of points
    % (note: this if-else structure is not particularly pretty)

    E = S;
    % read center, shift ellipsoid
    c = E.q; E = E + (-c);
    
    % compute rank and dimension
    r = rank(E); n = dim(E);
    
    % determine degeneracy: if so, project on proper subspace (via SVD)
    n_rem = n-r;
    [T,~,~] = svd(E.Q);
    E = T'*E;
    E = project(E,1:r);
    G = inv(sqrtm(E.Q));
    E = G*E;
end

% obtain center
c = center(E);

% quantile function for probability p of the chi-squared distribution
quantileFctValue = chi2inv(pr,dim(E));

% obtain covariance matrix
Sigma = E.Q/quantileFctValue;

% create N samples
x = zeros(dim(S),N);
remainingSamples = N;
idx = 1;

while remainingSamples > 0
    % create remaining number of samples of normal distribution
    pt = mvnrnd(c,Sigma,remainingSamples)';
    
    % check containment
    ptInside = contains(S,pt);
    nrInside = nnz(ptInside);
    remainingSamples = remainingSamples - nrInside;

    % store the ones that are contained, repeat loop for remaining number
    x(:,idx:idx+nrInside-1) = pt(:,ptInside');
    idx = idx+nrInside;
end

% degenerate ellipsoids: stack again, backtransform and shift
if isa(S,'ellipsoid')
    x = T*[inv(G)*pt;zeros(n_rem,N)] + c;
end

% ------------------------------ END OF CODE ------------------------------
