function I = interval(E)
% interval - Over-approximates an ellipsoid by an interval
%
% Syntax:
%    I = interval(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    I - interval object
%
% Example: 
%    E = ellipsoid([3 -1; -1 1],[1;0]);
%    I = interval(E);
%
%    figure; hold on;
%    plot(E);
%    plot(I,[1,2],'r');
%
% Other m-files required: interval (zonotope)
% Subfunctions: none
% MAT-files required: none
%
% See also: vertices, polytope

% Authors:       Victor Gassmann
% Written:       13-March-2019
% Last update:   04-July-2022 (VG, input checks)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check inputs
inputArgsCheck({{E,'att','ellipsoid','scalar'}});

n = dim(E);
E0 = ellipsoid(E.Q,zeros(size(E.q)));
dI = zeros(n,1);
Idty = eye(n);

% compute width of the ellipsoid in each dimension using support functions
for i=1:n
    dI(i) = supportFunc_(E0,Idty(:,i),'upper');
end

% construct the resulting interval
I = interval(-dI,dI) + E.q;

% ------------------------------ END OF CODE ------------------------------
