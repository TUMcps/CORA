function res = testLong_interval_rdivide
% testLong_interval_rdivide - unit_test_function of sine for
%    intervals, overloaded './' function for intervals
%
% Syntax:
%    res = testLong_interval_rdivide
%
% Inputs:
%    no
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Mark Wetzlinger
% Written:       08-August-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

numberRandTests = 5000;

% 1. random tests: division by scalar -------------------------------------

% division by zero returns -Inf,Inf
res_rand_scalar(1) = true;
for i=1:numberRandTests
    a = interval(-randi(1000,1,1)*rand(1),randi(1000,1,1)*rand(1));
    try
        b = a ./ 0; % error should be thrown here
        % in case no error is thrown, all should be -Inf|Inf
        if ~all(isinf(b.inf)) || ~all(isinf(b.sup))
            res_rand_scalar(1) = false;
            break
        end
    catch
        continue;
    end
end

% division by number other than zero should always produce a result
res_rand_scalar(2) = true;
for i=1:numberRandTests
    a = interval(-randi(1000,1,1)*rand(1),randi(1000,1,1)*rand(1));
    b = a ./ rand(1);
    if any(isnan(b.inf)) || any(isnan(b.sup)) || any(isinf(b.inf)) || any(isinf(b.sup))
        res_rand_scalar(2) = false;
        break
    end
end

% 2. random tests: division by interval -----------------------------------

% division by interval wholly negative or positive should always produce a result
res_rand_int(1) = true;
for i=1:numberRandTests
    randsign = sign(rand(1));
    randinf = randsign*randi(1000,1,1)*(1+rand(1));
    a = interval(randinf,randinf+rand(1));
    b = rand(1) ./ a;
    if any(isnan(b.inf)) || any(isnan(b.sup)) || any(isinf(b.inf)) || any(isinf(b.sup))
        res_rand_int(1) = false;
        break
    end
end

% division by interval containing 0 returns -Inf|Inf
res_rand_int(2) = true;
for i=1:numberRandTests
    a = interval(-randi(1000,1,1)*rand(1),randi(1000,1,1)*rand(1));
    try
        b = rand(1) ./ a; % should throw an error
        % in case no error is thrown...
        if ~all(isinf(b.inf)) || ~all(isinf(b.sup))
            res_rand_int(2) = false;
            break
        end
    catch
        continue
    end
end

% combine results
res = all(res_rand_scalar) && all(res_rand_int);

% ------------------------------ END OF CODE ------------------------------
