function x = randPoint(obj,varargin)
% randPoint - generates a random point within a given continuous set
%
% Syntax:  
%    x = randPoint(obj)
%    x = randPoint(obj,N)
%    x = randPoint(obj,N,type)
%    x = randPoint(obj,'all','extreme')
%
% Inputs:
%    obj - contSet object
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
% Last revision: ---

%------------- BEGIN CODE --------------

    % parse input arguments
    N = 1;
    type = 'gaussian';
    if nargin > 1 && ~isempty(varargin{1})
       if isa(varargin{1},'char')
           type = varargin{1};
       else
           N = varargin{1};
       end
    end
    if nargin > 2 && ~isempty(varargin{2}) && ~isa(varargin{1},'char')
       type = varargin{2}; 
    end

    % empty set
    if dim(obj) == 0
        [msg,id] = errEmptySet();
        error(id,msg);
    end
    
    % generate different types of extreme points according to 'type':
    % This function is only called directly when no randPoint
    % implementation for the subclass exists, or by a subclass itself
    % with type argument gaussian. Calling it directly with a type other than
    % gaussian therefore results in an error.
    if strcmp(type,'standard')
        error(strcat(errNoExactAlg(obj), " with randPoint type ''standard''"));
    elseif strcmp(type,'extreme')
        error(strcat(errNoExactAlg(obj), " with randPoint type ''extreme''"));
    elseif strcmp(type,'gaussian')
        
        % randPoint with type gaussian currently not supported for all
        % contSet subclasses, therefore check if subclass is supported
        if ~(isa(obj,'zonotope') || isa(obj,'interval') || isa(obj,'ellipsoid') || isa(obj,'mptPolytope'))
            error(strcat(errNoExactAlg(obj), " with randPoint type ''gaussian''"));
        end
        
        % get parameter pr
        if nargin > 2 && isa(varargin{2},'double')
            pr = varargin{2};
        elseif nargin > 3 && ~isempty(varargin{3})
            pr = varargin{3};
        else
            [msg,id] = errWrongInput('pr');
            error(id,msg);
        end
        
        % set is just a point
        if isa(obj,'zonotope') && isempty(generators(deleteZeros(obj)))
            x = center(obj); return;
        end
        
        % generates a random vector according to Gaussian distribution within a given set
        % enclose set by ellipsoid
        E = ellipsoid(obj);
        
        % obtain center
        c = center(E);
        
        % quantile function for probability p of the chi-squared distribution
        quantileFctValue = chi2inv(pr,dim(E));
        
        % obtain covariance matrix
        Sigma = E.Q/quantileFctValue;
        
        % create N samples
        x = zeros(dim(obj),N);
        
        for i = 1:N
            % create sample of normal distribution
            pt = mvnrnd(c,Sigma,1)';
            
            % correct value if not inclosed
            while ~(in(obj,pt))
                pt = mvnrnd(c,Sigma,1)';
            end
            x(:,i) = pt;
        end
    else
        [msg,id] = errWrongInput('type');
        error(id,msg);
    end    
%------------- END OF CODE --------------