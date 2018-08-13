function res = test_zonotope_guardProjection()
% test_zonotope_guardProjection - unit test function for computing the
% intersection of a reachable set with a halfspace when the initial set is
% a zonotope. This function is based on the HSCC'12 paper.
%
% Syntax:  
%    res = test_zonotope_guardProjection
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
% Written:      22-August-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% initial set
R0 = zonotope([0.7159, 0.0949, -0.0153, 0.0040, -0.0004; 0.9914, 0.0153, 0.0949, 0.0040, 0.0002]);

% guard
guard = halfspace([1;0],0.5);

% linear dynamics
A = [-1 -4; 4 -1];
B = [1;1];
twoDimSys=linearSys('twoDimSys',A,[1;1]); %initialize linear interval system

% location
reset.A=eye(2); %reset map for all transitions
reset.b=zeros(2,1); %reset map for all transitions
inv=interval([-10; -10], [10; 10]); %invariant for all locations
tran{1} = transition(guard,reset,2,'a','b'); 
loc = location('loc1',1,inv,tran,twoDimSys);   

% no partial intersection
partialIntersection = 0;

% options
options.timeStep = 0.04;
options.taylorTerms = 4; %number of taylor terms for reachable sets
options.uTrans = 0;
options.u = 0;
options.U = zonotope([0,0.1]); %input for reachability analysis
options.originContained = 1;
options.reductionTechnique = 'girard';
options.zonotopeOrder=10; %zonotope order
options.errorOrder=10;
options.x0 = center(R0);
options.maxProjectionError = inf;
options.linAlg = 1;

% call function
[R_guard, R_guard_noInt, splitFlag, RbeforeInt] = guardProjection(R0,guard,A,B,loc,partialIntersection,options);

% check if guard is close to a predefined interval
I = interval(R_guard{1});

%saved result
I_saved = interval( ...
        [0.5000000000000000; 0.8723179427185608], ...
        [0.5000000000000000; 1.2460156071355999]);
        
%check if slightly bloated versions enclose each other
res_1 = (I <= enlarge(I_saved,1+1e-8));
res_2 = (I_saved <= enlarge(I,1+1e-8));

%final result
res = res_1*res_2;


%------------- END OF CODE --------------
