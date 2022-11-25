function completed = example_nonlinear_reach_02_vanDerPol_polyZonotope()
% example_nonlinear_reach_02_vanDerPol_polyZonotope - example of non-linear 
%                                                     reachability analysis 
%
%    Example from [1] comparing reachability analysis for non-linear
%    systems using zonotopes and polynomial zonotopes for the Van-der-Pol
%    oscillator system.
%
% Syntax:  
%    completed = example_nonlinear_reach_02_vanDerPol_polyZonotope()
%
% Inputs:
%    no
%
% Outputs:
%    completed - boolean 
%
% References: 
%   [1] N. Kochdumper et al. "Sparse Polynomial Zonotopes: A Novel Set 
%       Representation for Reachability Analysis"

% Author:       Niklas Kochdumper
% Written:      02-January-2020
% Last update:  23-April-2020 (restucture params/options)
% Last revision:---


%------------- BEGIN CODE --------------

% Parameters --------------------------------------------------------------

params.tFinal = 6.74;                              % final time
R0 = zonotope([[1.4;2.4], diag([0.17, 0.06])]);    % initial set


% Reachability Settings ---------------------------------------------------

% settings
options.timeStep = 0.005;                           
options.taylorTerms = 4;                            
options.zonotopeOrder = 50;       
options.intermediateOrder = 50;
options.errorOrder = 20;

% reachability algorithm
options.alg = 'poly';
options.tensorOrder = 3;

% special settings for polynomial zonotopes
polyZono.maxDepGenOrder = 50;
polyZono.maxPolyZonoRatio = 0.01;
polyZono.restructureTechnique = 'reducePca';



% System Dynamics ---------------------------------------------------------

vanderPol = nonlinearSys(@vanderPolEq);



% Reachability Analysis (zonotope) ----------------------------------------
      
% adapted options
params.R0 = R0;

% compute reachable set
tic
R = reach(vanderPol, params, options);
tComp = toc;
disp(['computation time (zonotope): ',num2str(tComp)]);



% Reachability Analysis (polynomial zonotope) -----------------------------
      
% adapted options
params.R0 = polyZonotope(R0);

options.polyZono = polyZono;

% compute reachable set
tic
Rpoly = reach(vanderPol, params, options);
tComp = toc;
disp(['computation time (polynomial zonotope): ',num2str(tComp)]);



% Visualization -----------------------------------------------------------
    
figure; hold on;

% plot reachable set (zonotope)
handleLin = plot(R,[1,2],'FaceColor',[0 0.4 1],'EdgeColor','none');

% plot reachable set (polynomial zonotope)
handlePoly = plot(Rpoly,[1,2],'FaceColor',[0 0.8 0],'EdgeColor','none','Splits',2);

% plot initial set
plot(R0,[1,2],'w','Filled',true,'EdgeColor','k');

% label plot
xlabel('x_1');
ylabel('x_2');

legend([handleLin,handlePoly],'zonotope','sparse polynomial zonotope');
box on


% example completed
completed = 1;

%------------- END OF CODE --------------