function res = testLong_interval_cartProd
% testLong_interval_cartProd - unit test function of Cartesian
%    product for intervals
%
% Syntax:
%    res = testLong_interval_cartProd
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
% See also: none

% Authors:       Mark Wetzlinger
% Written:       03-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% number of tests
nrTests = 1000;

res = true;
for i=1:nrTests
    
    % random dimensions
    n1 = randi(10);
    n2 = randi(10);
    n3 = randi(10);

    % random lower and upper bounds
    lb1 = -rand(n1,1);
    ub1 = rand(n1,1);
    lb2 = -rand(n2,1);
    ub2 = rand(n2,1);
    num = randn(n3,1);

    % instantiate intervals
    I1 = interval(lb1,ub1);
    I2 = interval(lb2,ub2);

    % 1. interval-interval case

    % compute Cartesian product
    I_ = cartProd(I1,I2);

    % true result
    I = interval([lb1;lb2],[ub1;ub2]);

    % compare results
    if ~isequal(I,I_,1e-14)
        res = false; break
    end

    % 2. interval-numeric case

    % compute Cartesian product
    I_ = cartProd(I1,num);

    % true result
    I = interval([lb1;num],[ub1;num]);

    % compare results
    if ~isequal(I,I_,1e-14)
        res = false; break
    end

    % 2. interval-numeric case

    % compute Cartesian product
    I_ = cartProd(num,I1);

    % true result
    I = interval([num;lb1],[num;ub1]);

    % compare results
    if ~isequal(I,I_,1e-14)
        res = false; break
    end

end

% ------------------------------ END OF CODE ------------------------------
