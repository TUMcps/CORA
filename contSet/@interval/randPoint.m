function p = randPoint(obj,varargin)
% randPoint - computes random point in interval
%
% Syntax:  
%    p = randPoint(obj)
%    p = randPoint(obj,N)
%    p = randPoint(obj,type)
%    p = randPoint(obj,N,type)
%    p = randPoint(obj,'all','extreme')
%    p = randPoint(obj,N,'gaussian',pr)
%    p = randPoint(obj,'gaussian',pr)
%
% Inputs:
%    obj - interval object
%    N - number of random points
%    type - type of the random point ('standard', 'extreme', or 'gaussian')
%    pr - probability that a value is within the set (only type = 'gaussian')
%
% Outputs:
%    p - random point in interval
%
% Example: 
%    int = interval.generateRandom(2);
%    points = randPoint(int,20);
%
%    figure; hold on;
%    plot(int);
%    plot(points(1,:),points(2,:),'.k','MarkerSize',20);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/randPoint

% Author:       Mark Wetzlinger
% Written:      17-Sep-2019
% Last update:  25-June-2021 (MP, add type gaussian)
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    N = 1;
    type = 'standard';
    types = {'standard','extreme','gaussian'};
    
    if nargin > 1 && ~isempty(varargin{1})
        if (isnumeric(varargin{1}) && isscalar(varargin{1})) || ...
                (ischar(varargin{1}) && strcmp(varargin{1},'all'))
            % second input argument is number of points
            N = varargin{1};
            
            if nargin > 2
                if ischar(varargin{2}) && any(strcmp(varargin{2},types))
                    type = varargin{2};
                    if nargin > 3  && strcmp(type,'gaussian')
                        pr = varargin{3};
                    end
                end
            end
            
        elseif ischar(varargin{1}) && any(strcmp(varargin{1},types))
            % second input argument is type (sampling method)
            type = varargin{1};
            if nargin > 2 && strcmp(type,'gaussian')
                pr = varargin{2};
            end
            
        else
            [msg,id] = errWrongInput('N or type');
            error(id,msg);
        end
    end
    
    
    % generate different types of extreme points
    if strcmp(type,'gaussian')
        if nargin == 3
            p = randPoint@contSet(obj,type,pr);
        else
            p = randPoint@contSet(obj,N,type,pr);
        end
    else
        % get object properties
        c = center(obj); r = rad(obj); n = dim(obj);
        
        if isempty(obj)
            % empty set
            [msg,id] = errEmptySet(); error(id,msg);
        end
        
        if strcmp(type,'standard')
            
            p = c + (-1 + 2 * rand(length(r),N)) .* r;
            
            
        elseif strcmp(type,'extreme')
            
            % consider degenerate case
            ind = find(r > 0);
            if length(ind) < n
                obj = project(obj,ind);
                temp = randPoint(obj,N,type);
                p = c * ones(1,N);
                p(ind,:) = temp;
                return;
            end
            
            % return all extreme point
            if ischar(N) && strcmp(N,'all')
                
                p = vertices(obj);
                
                % generate random vertices
            elseif 10*N < 2^n
                
                p = zeros(n,N); cnt = 1;
                while cnt <= N
                    temp = sign(-1 + 2*rand(n,1));
                    if ~ismember(temp',p','rows')
                        p(:,cnt) = temp; cnt = cnt + 1;
                    end
                end
                p = c + p.*r;
                
                % select random vertices
            elseif N <= 2^n
                
                V = vertices(obj);
                ind = randperm(size(V,2));
                V = V(:,ind);
                p = V(:,1:N);
                
                % compute vertices and additional points on the boundary
            else
                
                V = vertices(obj);
                p = [V, zeros(n,N-size(V,2))];
                
                for i = size(V,2)+1:N
                    temp = sign(-1 + 2*rand(n,1));
                    ind = randi([1,n]);
                    temp(ind) = -1 + 2*rand();
                    p(:,i) = c + temp .* r;
                end
            end
        else
            [msg,id] = errWrongInput('type');
            error(id,msg);
        end
    end
end

%------------- END OF CODE --------------