function cPZ = conPolyZono(E)
% conPolyZono - Converts an ellipsoid to a constrained polynomial zonotope
%
% Syntax:
%    cPZ = conPolyZono(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    cPZ - conPolyZono object
%
% Example: 
%    E = ellipsoid([4 2;2 4],[1;1]);
%    cPZ = conPolyZono(E);
% 
%    figure; hold on; xlim([-2,4]); ylim([-2,4]);
%    plot(E,[1,2],'FaceColor','r');
%    plot(cPZ,[1,2],'b','Splits',20);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval

% Authors:       Niklas Kochdumper
% Written:       12-August-2019 
% Last update:   04-July-2022 (VG, avoid class array problems)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{E,'att','ellipsoid','scalar'}});

% eigenvalue decomposition of the ellipsoid matrix
[V,D] = eig(E.Q);

% dimension
n = dim(E);
% starting point
c = E.q;
% dependent generator matrix and exponent matrix
G = V*sqrt(D);
E = [eye(n);zeros(1,n)];

% constraints
A = [-0.5 ones(1,n)];
b = 0.5;
EC = [[zeros(n,1);1],[2*eye(n); zeros(1,n)]];
% identifiers
id = (1:n+1)';

% instantiate the constrained polynomial zonotope
cPZ = conPolyZono(c,G,E,A,b,EC,[],id);
    
% ------------------------------ END OF CODE ------------------------------
