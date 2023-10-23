function p = cubMapPoint(x1,x2,x3,T)
% cubMapPoint - evaluates the cubic map for points
%
% Syntax:
%    p = cubMapPoint(x1,x2,x3,T)
%
% Inputs:
%    x1,x2,x3 - nx1 vectors
%    T - nxn cell-array of tensors
%
% Outputs:
%    p - cubic map of points
%
% Example:
%    x1 = [1;-1]; x2 = [1;0]; x3 = [0;-1];
%    T{1,1} = [1 -1; 0 1]; T{1,2} = [1 0; 0 -1];
%    T{2,1} = [1 0; 0 -1]; T{2,2} = [-1 1; 0 -1];
% 
%    p = cubMapPoint(x1,x2,x3,T);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper
% Written:       16-August-2018
% Last update:   24-April-2023 (MW, moved here from unitTests)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get dimension for input argument check
n = size(x1,1);
% check input arguments
% inputArgsCheck({{x1,'att','numeric',{'nonnan','vector'}},...
%                 {x2,'att','numeric',{'nonnan','vector','size',[n,1]}},...
%                 {x3,'att','numeric',{'nonnan','vector','size',[n,1]}},...
%                 {T,'att','cell',{'size',[n,n]}}});

% init point
p = zeros(n,1);

% loop over all dimensions
for i = 1:n
    % loop over all quadratic matrices for this dimension
    for j = 1:size(T,2)
        p(i) = p(i) + (x1' * T{i,j} * x2) * x3(j);
    end
end

% ------------------------------ END OF CODE ------------------------------
