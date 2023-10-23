function p = quadMapPoint(x1,x2,Q)
% quadMapPoint - evaluates the quadratic map for points
%
% Syntax:
%    p = quadMapPoint(x1,x2,Q)
%
% Inputs:
%    x1,x2 - nx1 vectors
%    Q - nx1 cell-array of tensors
%
% Outputs:
%    p - quadratic map of points
%
% Example:
%    x1 = [1;-1]; x2 = [1;0];
%    Q{1,1} = [1 -1; 0 1]; Q{1,2} = [1 0; 0 -1];
% 
%    p = quadMapPoint(x1,x2,Q);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       24-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get dimension for input argument check
n = size(x1,1);
% check input arguments
% try
%     inputArgsCheck({{x1,'att','numeric',{'nonnan','vector'}},...
%                     {x2,'att','numeric',{'nonnan','vector','size',[n,1]}},...
%                     {Q,'att','cell',{'size',[n,1]}}});
% catch
%     inputArgsCheck({{x1,'att','numeric',{'nonnan','vector'}},...
%                     {x2,'att','numeric',{'nonnan','vector','size',[n,1]}},...
%                     {Q,'att','cell',{'size',[1,n]}}});
% end

% init point
p = zeros(n,1);
% loop over all dimensions
for i = 1:n
    p(i) = x1' * Q{i} * x2;
end

% ------------------------------ END OF CODE ------------------------------
