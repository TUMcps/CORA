function p = randPoint_(pgon, N, type, varargin)
% randPoint_ - draw random points within the polygon
%
% Syntax:
%    p = randPoint_(pgon, N, type)
%
% Inputs:
%    pgon - polygon
%    N - numeric, number of points
%    type - 'standard' or 'extreme'
%
% Outputs:
%    p - drawn points
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/randPoint

% Authors:       Niklas Kochdumper
% Written:       13-March-2020
% Last update:   11-October-2024 (TL, renamed to randPoint_, restructure)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% different types of random points
switch type
    case 'standard'
        p = aux_randPointStandard(pgon,N);

    case 'extreme'
        p = aux_randPointExtreme(pgon,N);

    otherwise
        throw(CORAerror('CORA:wrongValue','third',{'standard','extreme'}))
end
end


% Auxiliary functions -----------------------------------------------------

function p = aux_randPointStandard(pgon,N)
    % type = 'standard'

    p = zeros(2, N);
    
    list = triangulation(pgon);
    cnt = 1;
    cnt_ = 1;
    
    while cnt <= N
    
        % get random point by interpolation between the vertices
        V = vertices(list{cnt_});
        d = rand(3, 1);
        d = d ./ sum(d);
        p(:, cnt) = sum(V*diag(d), 2);
    
        % update counter
        cnt = cnt + 1;
        cnt_ = cnt_ + 1;
        if cnt_ > length(list)
            cnt_ = 1;
        end
    end
end

function p = aux_randPointExtreme(pgon,N)
    % type = 'extreme'

    % get all vertices
    V = vertices_(pgon);

    % filter nan values (multiple regions)
    V = V(:,~any(isnan(V),2));

    % return all extreme point
    if ischar(N) && strcmp(N, 'all')
        p = V;
        return
    end

    % check how many points should be returned
    numV = size(V,2);
    if N <= numV
        % take N vertices
        ind = randperm(numV);
        ind = ind(1:N);
        p = V(:,ind);
    else
        % sample additional vertices at boundary
        numAdd = N-numV;

        % select two subsequent points
        ind = randi(numV-1,[1,numAdd]);

        % compute convex combination
        d = rand(2,numAdd);
        d = d ./ sum(d);            
        Vadd = V(:,ind) .* d(1,:) +  V(:,ind+1) .* d(2,:);

        % output samples points
        p = [V,Vadd];
    end
    
end

% ------------------------------ END OF CODE ------------------------------
