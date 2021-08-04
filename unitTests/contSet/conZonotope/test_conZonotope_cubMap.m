function res = test_conZonotope_cubMap
% test_conZonotope_cubMap - unit test function for cubic multiplication of 
%                           constrained zonotopes
%
% Syntax:  
%    res = test_conZonotope_cubMap
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Niklas Kochdumper
% Written:      30-October-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = false;

%% ANALYTICAL TESTS

% TEST 1: Mixed Multiplication

% define zonotope
Z = [0 1 -1; 1 2 0];
A = [1 1];
b = 0;
cZ = conZonotope(Z,A,b);

% define third-order tensor
temp = [1 -1; 0 2];
T{1,1} = temp;
T{1,2} = temp;
T{2,1} = temp;
T{2,2} = temp;

% compute cubic map
cZres = cubMap(cZ,cZ,cZ,T);

% define ground truth
temp = [2 3 1 4 7 6 1 -1 -2 1 9 3 -3 12 -1 -4 21 3 -3 -7 3 -1 1 -1];
Z_ = [temp;temp];
A_ = zeros(3,23);
A_(1,1) = 1;
A_(1,2) = 1;
A_(2,3) = 1;
A_(3,5) = 1;
A_(3,8) = 1;
b_ = [0;0;0];

% check for correctness
if any(any(Z_-cZres.Z)) || any(any(A_-cZres.A)) || any(b_-cZres.b)
    error('conZonotope/cubMap: analytical test (mixed mul.) failed!');
end



% TEST 2: Cubic Multiplication

% define zonotope
Z = [0 1 -1; 1 2 0];
A = [1 1];
b = 0;
cZ = conZonotope(Z,A,b);

% define third-order tensor
temp = [1 -1; 0 2];
T{1,1} = temp;
T{1,2} = temp;
T{2,1} = temp;
T{2,2} = temp;

% compute cubic map
cZres = cubMap(cZ,T);

% define ground truth
temp = [16 13 -1 14 -4 21 -7 3 -1];
Z_ = [temp;temp];
A_ = [1 1 0 0 0 0 0 0];
b_ = 0;

% check for correctness
if any(any(Z_-cZres.Z)) || any(any(A_-cZres.A)) || any(b_-cZres.b)
    error('conZonotope/cubMap: analytical test failed!');
end

% no errors -> test successful
res = true;

end


% Auxiliary Functions -----------------------------------------------------

function p = cubMulPoint(x1,x2,x3,T)

    p = zeros(size(x1));
    
    % loop over all dimensions
    for i = 1:length(p)
       
         % loop over all quadratic matrices for this dimension
         for j = 1:size(T,2)
            p(i) = p(i) + (x1' * T{i,j} * x2) * x3(j);
         end
    end
end

%------------- END OF CODE --------------