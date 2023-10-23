function res = example_stl_bouncingBall
% example_stl_bouncingBall - example of signal temporal logic
% checking of the bouncing ball model
%
% Syntax:
%    res = example_stl_bouncingBall()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean

% Authors:       Benedikt Seidl
% Written:       04-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% Warning
disp('Warning: This example typically runs for many hours. Please comment the following lines to run it.');
completed = false;
return

alg = {'signals','sampledTime','rtl'}; %  % rtl takes a really long time

% Parameter ---------------------------------------------------------------

params.R0 = zonotope([1;0],diag([0.05,0.05]));
params.startLoc = 1;
params.tFinal = 2;


% Reachability Settings ---------------------------------------------------

% continuous dynamics
options.timeStep = 0.01;
options.taylorTerms = 10;
options.zonotopeOrder = 20;

% hybrid dynamics
options.guardIntersect = 'polytope';
options.enclose = {'box'};


% System Dynamics ---------------------------------------------------------

HA = bouncing_ball(-0.75);


% Reachability Analysis ---------------------------------------------------

tic
R = reach(HA, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


% Specification -----------------------------------------------------------

x = stl('x', 2);

phi = {
    globally(x(1) > 0.4, interval(0,0.1))
    finally(x(1) < 0.5, interval(0,0.5))
    finally(x(1) < 0.1, interval(0,1))
    globally(finally(x(2) < 0,interval(0,0.5)), interval(0,1.5))
    globally(finally(x(1) > 0.4, interval(0,1)), interval(0,0.5))
    globally(implies(x(1) < 0.1, finally(x(1) > 0.1, ...
        interval(0,0.5))), interval(0,1))
    finally(until(x(2) < 0, x(1) < 0.2, ...
        interval(0,1)), interval(0,0.1))
};


% Verification ------------------------------------------------------------

res = true;

for j = 1:length(alg)
    disp('-');
    disp(['checking alg ',alg{j}]);

    tFull = 0;

    for i = 1:length(phi)
        tic
        valid = modelChecking(R,phi{i},alg{j});
        tComp = toc;

        disp(['computation time for ',num2str(i),' is ',num2str(tComp)]);

        tFull = tFull + tComp;

        if ~valid
            disp('false negative');
        end

        res = res && valid;
    end

    disp(['computation time with ',alg{j},': ',num2str(tFull)]);
end

end

% ------------------------------ END OF CODE ------------------------------
