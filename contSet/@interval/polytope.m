function P = polytope(obj,varargin)
% polytope - Converts an interval object to a polytope object
%
% Syntax:  
%    P = polytope(ob,options)
%
% Inputs:
%    obj - interval object
%    options - options struct 
%
% Outputs:
%    P - polytope object
%
% Example: 
%    I = interval([1; -1], [2; 1]);
%    P = polytope(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope, interval

% Author:       Matthias Althoff
% Written:      22-July-2016 
% Last update:  14-June-2018
% Last revision:---

%------------- BEGIN CODE --------------

% generate polytope representation A*x <= b
leftLimit = obj.inf;
rightLimit = obj.sup;
dims = length(leftLimit);
A = zeros(0,dims);
b = zeros(0,1);

% restate each limit as an inequation a*x <= b_ (ignore unbounded limits)
for d = 1:dims
    if leftLimit(d) > -inf
        % form inequation "-x(d) <= -inf"
        a = zeros(1,dims);
        a(d) = -1;
        b_ = -1*leftLimit(d);
        
        % add row to polytope
        A = vertcat(A,a);
        b = vertcat(b,b_);
    end
    if rightLimit(d) < inf
        % form inequation "x(d) <= sup"
        a = zeros(1,dims);
        a(d) = 1;
        b_ = rightLimit(d);
        
        % add row to polytope
        A = vertcat(A,a);
        b = vertcat(b,b_);
    end
end

% call polytope constructor
param_struct.A = A;
param_struct.b = b;
P = mptPolytope(param_struct);

end
%------------- END OF CODE --------------