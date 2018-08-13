function res = commonPoint(obj, h_other)
% commonPoint - find arbitrary common point of two halfspaces
%
% Syntax:  
%    res = commonPoint(obj, h_other)
%
% Inputs:
%    obj - halfspace object
%    h_other - other halfspace object
%
% Outputs:
%    res - common point
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      27-August-2013
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%dimension
dim = length(obj.c);

%unit vector as initial point
x_0 = [1; zeros(dim-1,1)];

%first direction multiplier alpha_1
alpha_1 = obj.d - obj.c.'*x_0/(obj.c.'*obj.c);

%first projection
x_1 = x_0 + alpha_1*obj.c;

%new direction
n_2 = h_other.c - h_other.c.'*obj.c/norm(obj.c)^2*obj.c;

%second direction multiplier alpha_2
alpha_2 = h_other.d - h_other.c.'*x_1/(h_other.c.'*n_2);

%second projection
res = x_1 + alpha_2*n_2;



%------------- END OF CODE --------------