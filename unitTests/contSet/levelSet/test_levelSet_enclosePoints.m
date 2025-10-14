function res = test_levelSet_enclosePoints
% test_levelSet_enclosePoints - unit test function of enclosePoints
%
% Syntax:
%    res = test_levelSet_enclosePoints
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
% See also: enclosePoints

% Authors:       Niklas Kochdumper
% Written:       28-July-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 0.01;

% generate random point cloud
n = 1 + randi(4);
N = 10 + randi(20);
s = randi(10);

points = s*(-1 + 2*rand(n,N)) + s*(-1 + 2*rand(n,1));

% identify level set with different methods and polynomial orders
order = 2:4;
method = {'single','multiple'};

for i = 1:length(method)
    for j = 1:length(order)
        
        % identify enclosing level set
        ls = levelSet.enclosePoints(points,method{i},order(j));

        % check the results
        assert(all(contains(ls,points,'exact',tol)));
    end
end

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
