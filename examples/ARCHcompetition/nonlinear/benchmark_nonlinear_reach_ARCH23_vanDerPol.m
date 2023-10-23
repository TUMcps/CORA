function text = benchmark_nonlinear_reach_ARCH23_vanDerPol()
% benchmark_nonlinear_reach_ARCH23_vanDerPol - example of nonlinaer
%    reachability analysis. The guard intersections for the
%    pseudoinvariant are calculated with Girard's method
%
% Syntax:
%    benchmark_nonlinear_reach_ARCH23_vanDerPol()
%
% Inputs:
%    ---
%
% Outputs:
%    res - true/false
%

% Authors:       Mark Wetzlinger
% Written:       06-July-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

params.tFinal = 7;
params.startLoc = 1;

% original interval for x1,y1,x2,y2: interval([1.25;2.35],[1.55;2.45]);
% original interval for b: 1 to 3

% splits that work
% 1.25 to 1.4, -, 1 to 1.4
% 1.25 to 1.4, -, 1.4 to 1.8
% 1.25 to 1.4, -, 1.8 to 2.1
% 1.25 to 1.4, -, 2.1 to 2.3
% 1.25 to 1.4, -, 2.3 to 2.45
% 1.25 to 1.4, -, 2.45 to 2.66
% 1.25 to 1.4, -, 2.66 to 3

% 1.4 to 1.55, -, 1 to 1.4
% 1.4 to 1.55, -, 1.4 to 1.8
% 1.4 to 1.55, -, 1.8 to 2.2
% 1.4 to 1.55, -, 2.2 to 2.43
% 1.4 to 1.55, -, 2.43 to 2.66
% 1.4 to 1.55, -, 2.66 to 3

% split initial sets
I1 = interval([1.25;2.35],[1.4;2.45]);
I2 = interval([1.4;2.35],[1.55;2.45]);
B1 = interval(1,1.4);
B2 = interval(1.4,1.8);
B3 = interval(1.8,2.1);
B4 = interval(1.8,2.2);
B5 = interval(2.1,2.3);
B6 = interval(2.2,2.43);
B7 = interval(2.3,2.45);
B8 = interval(2.43,2.66);
B9 = interval(2.45,2.66);
B10 = interval(2.66,3);

R0{1} = cartProd(cartProd(I1,I1),B1);
R0{2} = cartProd(cartProd(I1,I1),B2);
R0{3} = cartProd(cartProd(I1,I1),B3);
R0{4} = cartProd(cartProd(I1,I1),B5);
R0{5} = cartProd(cartProd(I1,I1),B7);
R0{6} = cartProd(cartProd(I1,I1),B9);
R0{7} = cartProd(cartProd(I1,I1),B10);
R0{8} = cartProd(cartProd(I2,I2),B1);
R0{9} = cartProd(cartProd(I2,I2),B2);
R0{10} = cartProd(cartProd(I2,I2),B4);
R0{11} = cartProd(cartProd(I2,I2),B6);
R0{12} = cartProd(cartProd(I2,I2),B8);
R0{13} = cartProd(cartProd(I2,I2),B10);


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.005;
options.taylorTerms = 3;
options.zonotopeOrder = 100;

options.alg = 'poly';
options.tensorOrder = 3;
options.errorOrder = 10;
options.intermediateOrder = 50;

% polyZono.maxDepGenOrder = 50;
% polyZono.maxPolyZonoRatio = 0.02;
% polyZono.restructureTechnique = 'reducePca';

options.guardIntersect = 'zonoGirard';
options.enclose = {'pca'};
% options.verbose = true;


% System Dynamics ---------------------------------------------------------

sys = nonlinearSys(@coupledVanDerPol_ARCH23,5,1);


% Specifications ----------------------------------------------------------

% original specs: y1 <= 2.75, y2 <= 2.75 should be fulfilled at all times
speclim = 2.75;

hs1 = halfspace([0 -1 0 0 0],-speclim);
% spec = specification(hs1,'unsafeSet');
hs2 = halfspace([0 0 0 -1 0],-speclim);
spec = specification({hs1,hs2},'unsafeSet');


% Hybrid Automaton --------------------------------------------------------

% artificial guard set
c = [-0.9808 -0.7286 -0.9808 -0.7286 0]; d = -0.9811;
guard12 = conHyperplane(c,d);
inv1 = polytope(c,d);
guard23 = conHyperplane(-c,-d);
inv2 = polytope(-c,-d);
inv3 = polytope([1 0 0 0 0],1000);

% identity reset function
reset.A = eye(5); reset.c = zeros(5,1);

% instantiate transitions
trans12(1) = transition(guard12,reset,2);
trans23(1) = transition(guard23,reset,3);

% instantiate locations
loc(1) = location('loc1',inv1,trans12,sys);
loc(2) = location('loc2',inv2,trans23,sys);
loc(3) = location('loc3',inv3,transition(),sys);

% instantiate hybrid automaton for coupled van der Pol system
HA = hybridAutomaton(loc);


% Reachability Analysis ---------------------------------------------------

% simOptions.points = 20;
% simOptions.fracVert = 0.5;
% simRes = simulateRandom(sys,params,simOptions);

res = false(length(R0),1);
tic;
for i=1:length(R0)
    % display status
    disp("Split " + i + "/" + length(R0));

    % set initial set
    params.R0 = polyZonotope(R0{i});

    % reachability
    [R_,res(i)] = reach(HA,params,options,spec);

    % append reachSet objects
    if i == 1
        R = R_;
    else
        R = [R; R_];
    end
end
tComp = toc;
disp(['specifications verified: ', num2str(all(res))]);
disp(['computation time for van der Pol: ',num2str(tComp)]);


% Visualization -----------------------------------------------------------

% x1-y1 and x2-y2
figure;
projDim = {[1,2],[3,4]};

for i=1:length(projDim)
    subplot(1,2,i); hold on; box on;

    % axis limits
    xlim([-2.5,2.5]);
    ylim([-4.05,4.05]);
    
    % reachable set
    plot(R,projDim{i},'FaceColor',colorblind('b'));
    % last set
    plot(R(end).timePoint.set{end},projDim{i},'EdgeColor','k');
    % simulation
%     plot(simRes,projDim{i});
    
    % artificial guard sets
%     plot(guard,projDim{i},'k');

    % specs
    if i == 1
        plot(conHyperplane(hs1),projDim{i},'r--');
    elseif i == 2
        plot(conHyperplane(hs2),projDim{i},'r--');
    end
    
    % formatting
    xlabel(['x_' num2str(i)]);
    ylabel(['y_' num2str(i)]);
end

text = ['CVDP22, ,',num2str(all(res)),',',num2str(tComp),', , '];

% ------------------------------ END OF CODE ------------------------------
