function res = testLong_interval_supportFunc
% testLong_interval_supportFunc - unit test function of supportFunc
%
% Syntax:
%    res = testLong_interval_supportFunc
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       12-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;

% random cases
res = true;
nrOfTests = 1000;

for i=1:nrOfTests
    
    % random dimension
    n = randi(30);

    % init random full-dimensional interval
    lb = -rand(n,1);
    ub = rand(n,1);
    I = interval(lb,ub);

    % support function in random direction compared to zonotope
    dir = randn(n,1);
    [val,x] = supportFunc(I,dir,'upper');
    [valZ,xZ] = supportFunc(zonotope(I),dir,'upper');
    if ~withinTol(val,valZ,1e-14) || ~all(withinTol(x,xZ,1e-14))
        res = false; break
    end
    [val,x] = supportFunc(I,dir,'lower');
    [valZ,xZ] = supportFunc(zonotope(I),dir,'lower');
    if ~withinTol(val,valZ,1e-14) || ~all(withinTol(x,xZ,1e-14))
        res = false; break
    end

    % support function in axis aligned directions
    id = eye(n);
    for d=1:n
        dir = id(:,d);
        % upper direction
        s_upper = supportFunc(I,dir,'upper');
        % check with correct solution
        if abs(s_upper - I.sup(d)) > tol
            res = false; break;
        end

        % lower direction
        s_lower = supportFunc(I,dir,'lower');
        % check with correct solution
        if abs(s_lower - I.inf(d)) > tol
            res = false; break;
        end
    end

    % special cases...

    % support function of unit cube
    I = interval(-ones(n,1),ones(n,1));
    % simple direction (normalized!)
    dir = ones(n,1)/norm(ones(n,1));
    % support function
    vertex_upper = supportFunc(I,dir,'upper');
    vertex_lower = supportFunc(I,dir,'lower');
    % check with correct result
    if abs(vertex_upper - sqrt(n)) > tol || ...
            abs(vertex_lower + sqrt(n)) > tol
        res = false; break;
    end

    % interval is not full-dimensional
    randDim = randi(n);
    lb = -rand(n,1);
    ub = rand(n,1);
    lb(randDim) = 0;
    ub(randDim) = 0;
    I = interval(lb,ub);

    % support function in flattened dimension
    dir_zero = id(:,randDim);
    s_upper = supportFunc(I,dir_zero,'upper');
    s_lower = supportFunc(I,dir_zero,'lower');

    % check with correct result
    if s_upper ~= 0 || s_lower ~= 0
        res = false; break;
    end

end

% ------------------------------ END OF CODE ------------------------------
