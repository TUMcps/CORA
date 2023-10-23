function p = randPoint_(fs,N,type,varargin)
% randPoint_ - generates random points within a full-dimensional space
%    case R^0: only point is 0 (not representable in MATLAB)
%
% Syntax:
%    p = randPoint_(fs)
%    p = randPoint_(fs,N)
%    p = randPoint_(fs,N,type)
%    p = randPoint_(fs,'all','extreme')
%
% Inputs:
%    fs - fullspace object
%    N - number of random points
%    type - type of the random point ('standard', 'extreme')
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    O = fullspace(2);
%    p = randPoint(O);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/randPoint

% Authors:       Mark Wetzlinger
% Written:       05-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if fs.dimension == 0
    throw(CORAerror('CORA:notSupported','Sampling of R^0 not supported'));
end

if strcmp(type,'standard')
    % use built-in random sampling
    p = randn(fs.dimension,N);

elseif strcmp(type,'extreme')

    if ischar(N) && strcmp(N,'all')
        % sample all 2^n 'vertices', i.e., -Inf/+Inf
        % may throw an error if the dimension is too large...
        p = vertices(interval(fs));

    else
        % init with random points...
        p = randn(fs.dimension,N);
        % at least one entry has to be -Inf/+Inf
        pos_inf = ceil(rand(1,N)*fs.dimension);
        s = sign(randn(1,N));
        linear_idx = (0:fs.dimension:fs.dimension*(N-1)) + pos_inf;
        p(linear_idx) = s*Inf;

    end

end

% ------------------------------ END OF CODE ------------------------------
