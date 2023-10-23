function res = testLong_zonotope_intersectStrip
% testLong_zonotope_intersectStrip - unit test function of 
% intersectStrip
%
% Syntax:
%    res = testLong_zonotope_intersectStrip
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

% Authors:       Matthias Althoff
% Written:       08-September-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check whether the over-approximation of intersectStrip encloses the
% exact result

%% test methods handling multiple strips at a time
% specify strips
C = [1 0; 0 1; 1 1];
phi = [5; 3; 3];
y = [-2; 2; 2];


% polytope of intersected strips
P = polytope([1 0;-1 0; 0 1;0 -1; 1 1;-1 -1],[3;7;5;1;5;1]);

% specify zonotope
Z = zonotope([1 2 2 2 6 2 8; 1 2 2 0 5 0 6 ]);

% obtain over-approximative zonotope after intersection
Z_over{1} = intersectStrip(Z,C,phi,y,'normGen');
Z_over{2} = intersectStrip(Z,C,phi,y,'svd');
Z_over{3} = intersectStrip(Z,C,phi,y,'radius');
% flag 'alamo-volume' is disabled for non-scalar phi:
% Z_over{4} = intersectStrip(Z,C,phi,y,'alamo-volume');

% obtain exact solution
P_exact = polytope(Z) & P;

%% testing methods for single strips
% specify strip
C_single = [1 0];
phi_single = 5;
y_single = -2;

% polytope of intersected strips
P_single = polytope([1 0;-1 0],[3;7]);

% obtain over-approximative zonotope after intersection
Z_over_single{1} = intersectStrip(Z,C_single,phi_single,y_single,'normGen');
Z_over_single{2} = intersectStrip(Z,C_single,phi_single,y_single,'svd');
Z_over_single{3} = intersectStrip(Z,C_single,phi_single,y_single,'radius');
Z_over_single{4} = intersectStrip(Z,C_single,phi_single,y_single,'alamo-volume');
Z_over_single{5} = intersectStrip(Z,C_single,phi_single,y_single,'alamo-FRad');
Z_over_single{6} = intersectStrip(Z,C_single,phi_single,y_single,'bravo');

% obtain exact solution
P_exact_single = polytope(Z) & P_single;

%% exact solution enclosed in over-approximation?
% specify tolerance
tol = 1e-6;
Z_tol = zonotope([zeros(2,1),tol*eye(2)]);
% init resVec
resVec = zeros(length(Z_over)+length(Z_over_single),1);
% case for multiple strips
for i = 1:length(Z_over)
    resVec(i) = contains(polytope(Z_over{i} + Z_tol),P_exact);
end
% case for single strip
for i = 1:length(Z_over_single)
    resVec(i+length(Z_over)) = contains(polytope(Z_over_single{i} + Z_tol),P_exact_single);
end

% all methods over-aproximative?
res = all(resVec);
 
% figure; hold on 
% plot(Z,[1 2],'r-+');
% plot(P,[1 2],'r-*');
% plot(Z_over{1},[1 2],'b-+');
% plot(Z & P,[1 2],'g');
% plot(P_exact,[1 2],'b-*');
% legend('zonotope','strips','zonoStrips','zono&poly','exact');

% figure; hold on 
% plot(Z,[1 2],'r-+');
% plot(P_single,[1 2],'r-*');
% plot(Z_over_single{1},[1 2],'b-+');
% plot(Z_over_single{5},[1 2],'g-');
% plot(P_exact_single,[1 2],'b-*');
% legend('zonotope','strips','zonoStrips','bravoMethod','exact');


% ------------------------------ END OF CODE ------------------------------
