function res = testLong_conZonotope_conIntersect
% testLong_conZonotope_conIntersect - unit test function for the
%    constrained intersection of two constrained zonotope objects
%
% Syntax:
%    res = testLong_conZonotope_conIntersect
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Authors:       Matthias Althoff
% Written:       04-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% first constrained zonotope Z
Z = conZonotope(zeros(6,1),eye(6));

% second constrained zonotope Y
Y = conZonotope([-3.9389; -0.1408; -0.3043; 0.1978], 0.1*eye(4));

% constraint matrix R
R =  [...
    0    15     0     0     0     0; ...
    0     0     1     0     0     0; ...
    0     0     0     1     0     0; ...
    0     0     0     0     1     0];

% % sanity check 1
% intersection = (R*Z & Y);
% Z2 = R*Z;
% res = isIntersecting(Z2, Y);
% res = isIntersecting(R*Z, Y);

% intersection
Zint = conIntersect(Z, Y, R);
representsa_(Zint,'emptySet',eps);

% plot sets (for debugging)
% figure; 
% subplot(1,2,1); hold on
% plot(R*Z,[1 2],'b');
% plot(Y,[1 2],'r');
% subplot(1,2,2); hold on
% plot(R*Z,[3 4],'b');
% plot(Y,[3 4],'r');
% 
% figure;
% subplot(1,2,1); hold on
% plot(Z,[1 2],'b');
% plot(Zint,[1 2],'r');
% subplot(1,2,2); hold on
% plot(Z,[3 4],'b');
% plot(Zint,[3 4],'r');


% first constrained zonotope Z
Z_zono = zonotope([zeros(6,1),eye(6)]);

% second constrained zonotope Y
Y_zono = zonotope([[10.1944; -0.4028; 0.1440; 0.5506], 0.1*eye(4)]);

res2 = isIntersecting(R*Z_zono, Y_zono);

res3 = isIntersecting(polytope(R*Z_zono), polytope(Y_zono));

% the result below is an polytope (?)
res4 = polytope(R*Z_zono) & polytope(Y_zono);

% combine checks
res = res2 && res3;

% ------------------------------ END OF CODE ------------------------------
