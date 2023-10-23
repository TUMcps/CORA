function res = testLong_conZonotope_intersectStrip
% testLong_conZonotope_intersectStrip - unit test function of 
%    intersectStrip
%
% Syntax:
%    res = testLong_conZonotope_intersectStrip
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
% Written:       05-March-2021
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

% specify constrained zonotope cont here
Z = [1 2 2 2 6 2 8; 1 2 2 0 5 0 6];
A = [1 0 1 0 0 0];
b = 1;
cZ = conZonotope(Z,A,b);

% obtain over-approximative zonotope after intersection
cZ_over{1} = intersectStrip(cZ,C,phi,y,'normGen');
cZ_over{2} = intersectStrip(cZ,C,phi,y,'svd');
cZ_over{3} = intersectStrip(cZ,C,phi,y,'radius');

% obtain exact solution
P_exact = polytope(cZ) & P;

%% testing methods for single strips
% specify strip
C_single = [1 0];
phi_single = 5;
y_single = -2;

% polytope of intersected strips
P_single = polytope([1 0;-1 0],[3;7]);

% obtain over-approximative zonotope after intersection
cZ_over_single{1} = intersectStrip(cZ,C_single,phi_single,y_single,'normGen');
cZ_over_single{2} = intersectStrip(cZ,C_single,phi_single,y_single,'svd');
cZ_over_single{3} = intersectStrip(cZ,C_single,phi_single,y_single,'radius');
cZ_over_single{4} = intersectStrip(cZ,C_single,phi_single,y_single,'alamo-FRad');
cZ_over_single{5} = intersectStrip(cZ,C_single,phi_single,y_single,'bravo');

% obtain exact solution
P_exact_single = polytope(cZ) & P_single;

%% exact solution enclosed in over-approximation?
% specify tolerance
tol = 1e-6;
Z_tol = zonotope([zeros(2,1),tol*eye(2)]);
% init resVec
resVec = zeros(length(cZ_over)+length(cZ_over_single),1);
% case for multiple strips
for i = 1:length(cZ_over)
    resVec(i) = contains(polytope(cZ_over{i} + Z_tol),P_exact);
end
% case for single strip
for i = 1:length(cZ_over_single)
    resVec(i+length(cZ_over)) = ...
        contains(polytope(cZ_over_single{i} + Z_tol),P_exact_single,'exact',tol);
end

% all methods over-aproximative?
res = all(resVec);

% figure; hold on 
% plot(cZ,[1 2],'r-+');
% plot(P,[1 2],'r-*');
% plot(cZ_over{1},[1 2],'b-+');
% plot(cZ & P,[1 2],'g');
% plot(P_exact,[1 2],'b-*');
% legend('zonotope','strips','zonoStrips','zono&poly','exact');

% figure; hold on 
% plot(cZ,[1 2],'r-+');
% plot(P_single,[1 2],'r-*');
% plot(cZ_over_single{1},[1 2],'b-+');
% plot(cZ_over_single{5},[1 2],'g-');
% plot(P_exact_single,[1 2],'b-*');
% legend('zonotope','strips','zonoStrips','bravoMethod','exact');

% ------------------------------ END OF CODE ------------------------------
