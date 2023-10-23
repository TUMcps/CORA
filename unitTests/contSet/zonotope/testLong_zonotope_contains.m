function res = testLong_zonotope_contains
% testLong_zonotope_contains - unit test function of contains
%
% Syntax:
%    res = testLong_zonotope_contains
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

% Authors:       Matthias Althoff, Adrian Kulmburg
% Written:       26-July-2016
% Last update:   14-September-2019
%                01-July-2021 (AK, more tests, and merged with testLong_zonotope_containsPoint)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% Basic Containment Test
% create zonotope
Z1 = zonotope([-4, -3, -2, -1; 1, 2, 3, 4]);

% create parallelotopes
P1 = zonotope([-3.8, -4, 3; 1.2, 3, -4]);
P2 = zonotope([-3.8, -8, 2; 1.2, 10, -10]);

% obtain results
int_1 = contains(P1,Z1);
int_2 = contains(P2,Z1);

% true results
true_int_1 = 0;
true_int_2 = 1;   

% check result
res = (int_1==true_int_1) & (int_2==true_int_2);

%% Advanced Containment Test with all methods
in_venum_1 = contains(P1,Z1,'venum');
in_venum_2 = contains(P2,Z1,'venum');
res_venum = (in_venum_1==true_int_1) & (in_venum_2==true_int_2);

in_polymax_1 = contains(P1,Z1,'polymax');
in_polymax_2 = contains(P2,Z1,'polymax');
res_polymax = (in_polymax_1==true_int_1) & (in_polymax_2==true_int_2);

in_opt_1 = contains(P1,Z1,'opt',0,200);
in_opt_2 = contains(P2,Z1,'opt',0,200);
res_opt = (in_opt_1==true_int_1) & (in_opt_2==true_int_2);

in_st_1 = contains(P1,Z1,'st');
in_st_2 = contains(P2,Z1,'st');
res_st = (in_st_1==true_int_1) & (in_st_2==true_int_2);

res = res & res_venum & res_polymax & res_opt & res_st;

%% Point Containment
% create a zonotope
n = 2;
gens = 5;
% rand gives value of [0,1]
Z = zonotope([zeros(n,1),rand(n,gens)]);

% Single points -----------------------------------------------------------
% create a point inside: center of zonotope
p_inside = center(Z);
% ...and outside: add all generators (with max of rand -> 1)
p_outside = gens*ones(n,1);

% check if correct results for containment
res_inside  = contains(Z, p_inside);
res_outside = contains(Z, p_outside);

% point on the border should also count as outside (?)
% p_border = sum(generators(Z),2);
% res_border = containsPoint(Z, p_border);

res_single = res_inside && ~res_outside;% && ~res_border;
% -------------------------------------------------------------------------

%% Array of points
num = 10;
p_array = gens*(ones(n,num)+rand(n,num));
res_array = all(contains(Z,p_array));
% -------------------------------------------------------------------------

% add results
res = res && res_single && ~res_array;

% ------------------------------ END OF CODE ------------------------------
