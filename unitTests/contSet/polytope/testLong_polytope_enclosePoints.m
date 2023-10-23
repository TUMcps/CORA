function res = testLong_polytope_enclosePoints
% testLong_polytope_enclosePoints - unit test function of
%    enclosePoints
%
% Syntax:
%    res = testLong_polytope_enclosePoints
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       04-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 25;

for i=1:nrTests
    
    % random dimension
    n = randi(10);

    % generate random point cloud
    A = rand(n,n);
    
    c = randn(n,1);
    sigma = 0.5*(A+A') + n*eye(n);
    points = mvnrnd(c,sigma,2*n)';
    
    % generate random polytope
    P = polytope.enclosePoints(points);

    % check if polytope contains points
    if ~all(contains(P,points))
        res = false; return
    end

end

% ------------------------------ END OF CODE ------------------------------
