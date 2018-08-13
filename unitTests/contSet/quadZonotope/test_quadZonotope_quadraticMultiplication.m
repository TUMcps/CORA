function res = test_quadZonotope_quadraticMultiplication
% test_quadZonotope_quadraticMultiplication - unit test function of quadraticMultiplication
%
% Syntax:  
%    res = test_quadZonotope_quadraticMultiplication
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
% Written:      10-August-2017
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% create zonotope
c = [1;2];
G = [2 1;3 1];
Gquad = [1 2;3 5];
Gsquare = [3;4];
Grest = [4 5 1;2 3 1];
qZ1 = quadZonotope(c,G,Gquad,Gsquare,Grest);

% create matrices
Q{1} = [1 -2; -3 4];
Q{2} = [0.5 0; 2 -1];

% obtain result
qZ2 = quadraticMultiplication(qZ1,Q);
intComp = interval(qZ2);

% draw random points from zontope and compute approximately compute the
% quadratic mapping
N = 10000;

for j = 1:2
    value = zeros(N,1);
    for i = 1:N
        temp = randPoint(qZ1);
        value(i) = temp' * Q{j} * temp;
    end
    intReal(j,1) = interval(min(value),max(value));
end


% display the results
intReal
intComp

% check for correctness
res = 1;
for i = 1:2
   if supremum(intReal(i)) > supremum(intComp(i))
       res = 0;
       break;
   elseif infimum(intReal(i)) < infimum(intComp(i))
       res = 0;
       break;
   end
end

if res
    disp('test_zonotope_quadraticMultiplication successful');
else
    disp('test_zonotope_quadraticMultiplication failed');
end

%------------- END OF CODE --------------