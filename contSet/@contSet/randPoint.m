function x = randPoint(S,N,type,pr)
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
%    type - type of the random point ('standard', 'extreme', or 'gaussian')
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

% Author:        Matthias Althoff
% Written:       19-November-2020
% Last update:   19-June-2021 (MP, generalization)
%                08-March-2023 (MW, reduce number of contains-calls)
% Last revision: ---

%------------- BEGIN CODE --------------
   
% randPoint with type gaussian currently not supported for all contSet
% subclasses, therefore check if subclass is supported
if ~(isa(S,'zonotope') || isa(S,'interval') || isa(S,'ellipsoid') || isa(S,'mptPolytope'))
    throw(CORAerror('CORA:notSupported',"The function randPoint for " + class(S) + ...
        " does not support the type = 'gaussian'."));
end

% zonotope: set is just a point
if isa(S,'zonotope') && isempty(generators(deleteZeros(S)))
    x = center(S); return;
end

% generates a random vector according to Gaussian distribution within a given set
% enclose set by ellipsoid
E = ellipsoid(S);

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

%------------- END OF CODE --------------