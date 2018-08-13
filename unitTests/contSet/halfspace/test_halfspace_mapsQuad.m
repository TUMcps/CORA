function res = test_halfspace_mapsQuad()
% test_halfspace_mapsQuad - unit test function for computing the maps
% requires for a quadratic abstraction of the map from a reachable set to
% the one that intersects with a hyperplane as guard set. This function is 
% based on the HSCC'12 paper. We will also show the relationship between
% the used symbols in the function and the paper.
%
% Syntax:  
%    res = test_halfspace_mapsQuad
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

% Author:       Matthias Althoff
% Written:      23-August-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% halfspace
n = [2; -1]; % normal vector
d = -4; % distance
obj = halfspace(n,d);

% initial set
G = [1, -1; 0, 2];
R0 = zonotope([[4; 4], G]);

% initial state
x0 = center(R0);

% flow 
f = [-4; 2];

% system matrix
A = [-1 -4; 4 -1];

% call function
[k,L,Q,lambdaQuad_inv_zono,Lambda_int] = mapsQuad(obj, A, f, x0, R0);

% identity matrix
I = eye(length(n));

% values from the HSCC'12 paper for this example---------------------------
% auxiliary values
Lambda = -62;
Upsilon = [-6 -7];
Theta = [76; -118];
Omega = [0, 20; -20, 0];

% k vector (y_h after eq. (15) in HSCC'12 without x_0)
k_calc = 4/31*[-24; 14];

% L matrix (Proposition 2 in HSCC'12)
L_calc = 4/31*[-1 -4; 4 -1] + 1/961*[-456 708; 266 -413]; % identity removed from L compared to HSCC'12 for computational reasons

% Lambda_interval
Lambda_int_calc = interval(-76,-48);

% 1/Lambda_int_calc^2     
lambdaQuad = 1/Lambda_int^2;

% Theta_interval
Theta_int = interval([36; -158],[116; -78]); % depedndencies of y have been used

% Psi_c
ThetaDivLambda = interval([-29/12; 39/38],[-9/19; 79/24]); %IMPROVEMENT POSSIBLE HERE SINCE ONE CAN DIVIDE BETTER
Psi_c = [-659/456; 1969/912];
Psi_g = [443/456; 1033/912];

% Q tensor (Proposition 3 in HSCC'12)
% auxiliary values
Theta_c = [76; -118]; %=Theta
Theta_g{1} = [0; -20]; %(alpha = 1)
Theta_g{2} = [40; 20]; %(alpha = 2)

% center matrices
C1 = [5020/19, -3428/19; -6228/19, 4153/19];
C2 = [6939/19, -62843/114; -17681/38, 150289/228];

% first generator matrices
G1{1} = [-5316/19, -6202/19; -6198/19, -7231/19];
G2{1} = [3101/19, 21707/114; 7231/38, 50617/228];

% second round of generator matrices
G1{2} = [659/38, 4613/228; 1071/76, 59177/456];
G2{2} = [-1318/19, -4613/57; -1071/19, 18343/114];

G1{3} = [1573/38, -40669/228; -16823/76, -169441/456];
G2{3} = [8057/19, 1573/38; 173/38, -16823/76];

% final round of generator matrices

G1{4} = [-443/38, -3101/228; -1033/76, -7231/456];
G2{4} = [886/19, 3101/57; 1033/19, 7231/114];

G1{5} = [-3101/38, -21707/228; -7231/76, -50617/456];
G2{5} = [-1329/19, -3101/38; -3099/38, -7231/76];
%--------------------------------------------------------------------------

% compare results----------------------------------------------------------

% k vector
k_comp = all(abs(k - k_calc) < 1e-15);

% L matrix
L_comp = all(all((abs(L - L_calc) < 1e-15)));

% Q matrices
Q_comp(1) = all(all((abs(Q{1}.center - C1) < 1e-12)));
Q_comp(2) = all(all((abs(Q{2}.center - C2) < 1e-12)));
for i=1:length(G1)
    Q_comp(2+i) = all(all((abs(Q{1}.generator{i} - G1{i}) < 1e-12)));
end
for i=1:length(G2)
    Q_comp(2+length(G1)+i) = all(all((abs(Q{2}.generator{i} - G2{i}) < 1e-12)));
end

% Lambda interval
Lambda_comp = (Lambda_int == Lambda_int_calc);

% lambdaQuad_inv_zono
lambdaQuad_comp(1) = (abs(lambdaQuad_inv_zono.center - mid(lambdaQuad) < 1e-15));
lambdaQuad_comp(2) = (abs(lambdaQuad_inv_zono.generator{1} - rad(lambdaQuad) < 1e-15));
%--------------------------------------------------------------------------

res = k_comp*L_comp*all(Q_comp)*all(lambdaQuad_comp)*Lambda_comp;

%------------- END OF CODE --------------
