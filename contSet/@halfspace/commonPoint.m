function p = commonPoint(hs1,hs2)
% commonPoint - find arbitrary common point of two halfspaces
%
% Syntax:
%    p = commonPoint(hs1,hs2)
%
% Inputs:
%    hs1 - halfspace object
%    hs2 - halfspace object
%
% Outputs:
%    res - point that is an element of both halfspaces
%
% Example:
%    hs1 = halfspace([1 1],2);
%    hs2 = halfspace([1 -1],2);    
%    p = commonPoint(hs1,hs2);
% 
%    figure; hold on; xlim([-3,3]); ylim([-3,3]);
%    plot(hs1);
%    plot(hs2);
%    scatter(p(1),p(2),16,'r','filled');
% 
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       27-August-2013
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% dimension
n = dim(hs1);

% unit vector as initial point
x_0 = [1; zeros(n-1,1)];

% first direction multiplier alpha_1
alpha_1 = hs1.d - hs1.c.'*x_0/(hs1.c.'*hs1.c);

% first projection
x_1 = x_0 + alpha_1*hs1.c;

% new direction
n_2 = hs2.c - hs2.c.'*hs1.c/norm(hs1.c)^2*hs1.c;

% second direction multiplier alpha_2
alpha_2 = hs2.d - hs2.c.'*x_1/(hs2.c.'*n_2);

% second projection
p = x_1 + alpha_2*n_2;

% ------------------------------ END OF CODE ------------------------------
