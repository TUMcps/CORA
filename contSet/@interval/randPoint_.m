function p = randPoint_(I,N,type,varargin)
% randPoint_ - computes random point in interval
%
% Syntax:
%    p = randPoint_(I)
%    p = randPoint_(I,N)
%    p = randPoint_(I,N,type)
%    p = randPoint_(I,'all','extreme')
%
% Inputs:
%    I - interval object
%    N - number of random points
%    type - type of the random point ('standard', 'extreme')
%
% Outputs:
%    p - random point in interval
%
% Example: 
%    I = interval([-2;1],[3;2]);
%    p = randPoint(I,20);
%
%    figure; hold on;
%    plot(I);
%    plot(p(1,:),p(2,:),'.k','MarkerSize',10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/randPoint, zonotope/randPoint_

% Authors:       Mark Wetzlinger
% Written:       17-September-2019
% Last update:   25-June-2021 (MP, add type gaussian)
%                19-August-2022 (MW, integrate standardized pre-processing)
% Last revision: 27-March-2023 (MW, rename randPoint_)

% ------------------------------ BEGIN CODE -------------------------------

% get object properties
c = center(I); r = rad(I); n = dim(I);

% generate different types of extreme points
if strcmp(type,'standard') || startsWith(type,'uniform')
    
    if size(r,2) > 1
        if size(r,1) > 1
            % both dimensions larger than 1 -> interval matrix
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'interval/randPoint not defined for interval matrices!'));
        else
            % row interval
            p = c + (-1 + 2 * rand(N,length(r))) .* r;
        end
    else
        % column interval
        p = c + (-1 + 2 * rand(length(r),N)) .* r;
    end
    
elseif strcmp(type,'extreme')
    
    % consider degenerate case
    ind = find(r > 0);
    if length(ind) < n
        I = project(I,ind);
        temp = randPoint_(I,N,type);
        p = c * ones(1,N);
        p(ind,:) = temp;
        return;
    end
    
    % return all extreme points
    if ischar(N) && strcmp(N,'all')
        p = vertices(I);
        
    elseif 10*N < 2^n
        % generate random vertices
        p = zeros(n,N); cnt = 1;
        while cnt <= N
            temp = sign(-1 + 2*rand(n,1));
            if ~ismember(temp',p','rows')
                p(:,cnt) = temp; cnt = cnt + 1;
            end
        end
        p = c + p.*r;
        
        
    elseif N <= 2^n
        % select random vertices
        V = vertices(I);
        ind = randperm(size(V,2));
        V = V(:,ind);
        p = V(:,1:N);
        
    else
        % compute vertices and additional points on the boundary
        V = vertices(I);
        p = [V, zeros(n,N-size(V,2))];
        
        for i = size(V,2)+1:N
            temp = sign(-1 + 2*rand(n,1));
            ind = randi([1,n]);
            temp(ind) = -1 + 2*rand();
            p(:,i) = c + temp .* r;
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
