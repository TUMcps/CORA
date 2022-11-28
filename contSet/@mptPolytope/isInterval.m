function [res,I] = isInterval(P)
% isInterval - checks if a polytope can be equivalently represented by an 
%              interval
%
% Syntax:  
%    [res,I] = isInterval(P)
%
% Inputs:
%    P - mptPolytope object
%
% Outputs:
%    res - true/false
%    I - equivalent interval object
%
% Example: 
%    A = [1 0;-1 0;0 1;0 -1];
%    b = [3;-2;3;-2];
%    P = mptPolytope(A,b);
%
%    [res,I] = isInterval(P);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval

% Author:       Niklas Kochdumper
% Written:      26-November-2021 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = false; I = [];

% fast initial check
if ~all(sum(abs(sign(P.P.A)),2) == 1)
    return;
end

% construct equivalent interval
n = dim(P); l = -inf*ones(n,1); u = inf*ones(n,1);

for i = 1:n
   ind = find(P.P.A(:,i) ~= 0);
   for j = 1:length(ind)
       if P.P.A(ind(j),i) > 0
          u(i) = min(u(i),P.P.b(ind(j))/P.P.A(ind(j),i)); 
       else
          l(i) = max(l(i),P.P.b(ind(j))/P.P.A(ind(j),i));
       end
   end
end

% check if the interval is bounded
if ~any(isinf(u)) && ~any(isinf(l))
    I = interval(l,u); res = true; 
end

%------------- END OF CODE --------------