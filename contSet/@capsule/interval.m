function I = interval(C)
% interval - Over-approximate a capsule by an interval
%
% Syntax:
%    I = interval(C)
%
% Inputs:
%    C - capsule object
%
% Outputs:
%    I - interval object
%
% Example: 
%    C = capsule([0;0],[-2;2],2);
%    I = interval(C);
%
%    figure; hold on;
%    plot(C,[1,2],'k');
%    plot(I,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope

% Authors:       Niklas Kochdumper
% Written:       20-November-2019 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% initialization
n = length(C.c);
ub = zeros(n,1);
lb = zeros(n,1);
E = eye(n);

% loop over all dimensions
for i = 1:n
    ub(i) = supportFunc_(C,E(:,i),'upper');
    lb(i) = supportFunc_(C,E(:,i),'lower');        
end

% instantiate interval
I = interval(lb,ub);

% ------------------------------ END OF CODE ------------------------------
