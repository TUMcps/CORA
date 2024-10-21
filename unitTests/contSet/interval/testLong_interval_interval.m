function res = testLong_interval_interval
% testLong_interval_interval - unit test function of interval
%
% Syntax:
%    res = testLong_interval_interval
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
% Written:       20-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-12;

res = true;
nrOfTests = 1000;
for i=1:nrOfTests

    % random dimension
    n = randi(25);
    
    % random lower bound, random upper bound
    a = -rand(n,1);
    b = rand(n,1);
    a_mat = -rand(n+1);
    b_mat = rand(n+1);
    
    % admissible initializations
    I = interval(a,b);
    assertLoop(all(abs(I.inf - a) <= tol),i)
    assertLoop(all(abs(I.sup - b) <= tol),i)
    
    I = interval(a);
    assertLoop(all(abs(I.inf - a) <= tol),i)
    assertLoop(all(abs(I.sup - a) <= tol),i)
    
    I = interval(a_mat);
    assertLoop(all(abs(I.inf - a_mat) <= tol),i)
    assertLoop(all(abs(I.sup - a_mat) <= tol),i)
    
    I = interval(a_mat,b_mat);
    assertLoop(all(abs(I.inf - a_mat) <= tol),i)
    assertLoop(all(abs(I.sup - b_mat) <= tol),i)
    
    % wrong initializations
    a_large = 1+rand(n,1);
    b_small = -1-rand(n,1);
    a_plus1 = -rand(n+1,1);
    b_plus1 = rand(n+1,1);
    
    % lower limit larger than upper limit
    assertThrowsAs(@interval,'CORA:wrongInputInConstructor',a,b_small);
    assertThrowsAs(@interval,'CORA:wrongInputInConstructor',a_large,b);

    % size of limits do not match
    assertThrowsAs(@interval,'CORA:wrongInputInConstructor',a_plus1,b);
    assertThrowsAs(@interval,'CORA:wrongInputInConstructor',a,b_plus1);
    assertThrowsAs(@interval,'CORA:wrongInputInConstructor',a_mat,b);
    assertThrowsAs(@interval,'CORA:wrongInputInConstructor',a,b_mat);
    
    % too many input arguments
    assertThrowsAs(@interval,'CORA:numInputArgsConstructor',a,b,b);
end

% ------------------------------ END OF CODE ------------------------------
