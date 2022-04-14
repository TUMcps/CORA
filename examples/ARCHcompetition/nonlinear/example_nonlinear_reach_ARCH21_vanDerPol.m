function example_nonlinear_reach_ARCH21_vanDerPol()
% example_nonlinear_reach_ARCH21_vanDerPol - example of
% nonlinear reachability analysis. The guard intersections for the
% pseudoinvariant are calculated with Girards method
%
% Syntax:  
%    example_nonlinear_reach_ARCH21_vanDerPol()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% 
% Author:       Matthias Althoff
% Written:      02-April-2017
% Last update:  13-March-2018
% Last revision:---


%------------- BEGIN CODE --------------


% Parameter  --------------------------------------------------------------

% case mu = 1
params1.tFinal = 7;
params1.startLoc = 1;

int1 = interval([1.25;2.35],[1.55;2.45]);
params1.R0 = zonotope(cartProd(int1,int1));


% case mu = 2
params2.tFinal = 8;
params2.startLoc = 1;

int2 = interval([1.55;2.35],[1.85;2.45]);
params2.R0 = zonotope(cartProd(int2,int2));


% Reachability Settings ---------------------------------------------------

% case mu = 1
options1.taylorTerms = 3;
options1.zonotopeOrder = 20;

options1.alg = 'lin';
options1.tensorOrder=3;
options1.errorOrder = 10;
options1.intermediateOrder = 20;

options1.guardIntersect = 'zonoGirard';
options1.enclose = {'pca'};

options1.timeStep = 0.01; 

% case mu = 2
options2.taylorTerms = 3;
options2.zonotopeOrder = 100;

options2.alg = 'lin';
options2.tensorOrder=3;
options2.errorOrder = 10;
options2.intermediateOrder = 50;

options2.guardIntersect = 'zonoGirard';
options2.enclose = {'pca'};

options2.timeStep = 0.001; 


% System Dynamics ---------------------------------------------------------

sys1 = nonlinearSys(@coupledVanDerPol1);
sys2 = nonlinearSys(@coupledVanDerPol2);


% Specifications ----------------------------------------------------------

% case mu = 1
hs1 = halfspace([0 -1 0 0],-2.75);
hs2 = halfspace([0 0 0 -1],-2.75);

spec1 = specification({hs1,hs2},'unsafeSet');

% case mu = 1
hs1 = halfspace([0 -1 0 0],-4.05);
hs2 = halfspace([0 0 0 -1],-4.05);

spec2 = specification({hs1,hs2},'unsafeSet');


% Hybrid Automaton --------------------------------------------------------

% case mu = 1
c = [-0.9808 -0.7286 -0.9808 -0.7286]; d = -0.9811;
guard = conHyperplane(c,d);
inv1 = mptPolytope(c,d);
inv2 = mptPolytope([1 0 0 0],1000);

reset.A = eye(4); reset.b = zeros(4,1);

trans{1} = transition(guard,reset,2);

loc{1} = location('loc1',inv1,trans,sys1);
loc{2} = location('loc2',inv2,[],sys1);

HA1 = hybridAutomaton(loc);

% case mu = 2
c = [-0.4105 -0.1543 -0.4105 -0.1543]; d = -1.3;
guard = conHyperplane(c,d);
inv = mptPolytope(c,d);
reset.A = eye(4); reset.b = zeros(4,1);
trans{1} = transition(guard,reset,2);
loc{1} = location('loc1',inv,trans,sys2);

c = [0.4409 0.1886 0.4409 0.1886]; d = -1.2617;
guard = conHyperplane(c,d,[0 -1 0 0],0);
inv = mptPolytope([0 1 0 0],0.6);
reset.A = eye(4); reset.b = zeros(4,1);
trans{1} = transition(guard,reset,3);
loc{2} = location('loc2',inv,trans,sys2);

inv = mptPolytope([1 0 0 0],1000);
loc{3} = location('loc3',inv,[],sys2);

HA2 = hybridAutomaton(loc);


% Reachability Analysis ---------------------------------------------------

% reachability analysis for mu = 1
tic
[R1,res1] = reach(HA1,params1,options1,spec1);
tComp = toc;
disp(['specifications verified (\mu = 1): ', num2str(res1)]);
disp(['computation time for van der Pol (\mu = 1): ',num2str(tComp)]);

% reachability analysis for mu = 2
tic
[R2,res2] = reach(HA2,params2,options2,spec2);
tComp = toc;
disp(['specifications verified (\mu = 2): ', num2str(res2)]);
disp(['computation time for van der Pol (\mu = 2): ',num2str(tComp)]);


% Visualization -----------------------------------------------------------

figure; hold on; box on;

% case mu = 1
plot(R1,[1,2],'r','Filled',true,'EdgeColor','none');

% case mu = 2
plot(R2,[1,2],'b','Filled',true,'EdgeColor','none');

% formatting
xlim([-2.5,2.5]);
ylim([-4.05,4.05]);
xlabel('x_1');
ylabel('y_1');

%------------- END OF CODE --------------