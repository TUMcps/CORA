function completed = example_linear_reach_ARCH18_iss_Uconst_safe()
% example_linear_reach_ARCH18_iss_Uconst_safe - example for linear 
%       dynamics with constant inputs. It is checked if the safe
%       specification is satisfied
%
% Syntax:  
%    example_linear_reach_ARCH18_iss_Uconst_safe()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% 
% Author:       Niklas Kochdumpemr
% Written:      20-June-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% load system matrices
load('iss.mat');

% construct extended system matrices (inputs as additional states)
dim = length(A);
A_ = [A,B;zeros(size(B,2),dim + size(B,2))];
B_ = zeros(dim+size(B,2),1);
C_ = [C,zeros(size(C,1),size(B,2))];

% construct the linear system object
sys = linearSys('iss',A_,B_,[],C_);



% Options -----------------------------------------------------------------

R0 = [interval(-0.0001*ones(dim,1),0.0001*ones(dim,1));interval([0;0.8;0.9],[0.1;1;1])];
options.x0=mid(R0); % initial state for simulation

options.taylorTerms=10; % number of taylor terms for reachable sets
options.zonotopeOrder=10; % zonotope order
options.originContained=1;
options.reductionTechnique='girard';
options.linAlg = 1;
options.compOutputSet = 1;

options.verifySpecs = @specificationSafe;

options.uTrans=0; % center of input set
options.U=zonotope(0); % input for reachability analysis


% Reachability analysis ---------------------------------------------------

% compute reachable set using zonotopes
tic
options.R0=zonotope(R0); % initial state for reachability analysis
options.tStart=0; % start time
options.tFinal=20; % final time
options.timeStep = 0.02;
[~,~,~,Rout,~,res] = reachIss(sys,options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);

if res
   disp('Specification satisfied!'); 
else
   disp('Specification violated!'); 
end



% Simulation --------------------------------------------------------------

% create random simulations. RRTs would provide better results, but are
% computationally more demanding
runs = 10;
fractionVertices = 0.5;
fractionInputVertices = 0.5;
inputChanges = 20;
simRes = simulate_random(sys, options, runs, fractionVertices, fractionInputVertices, inputChanges);




% Visualization -----------------------------------------------------------

% plot the solution
hold on
i = 1;
iDim = 3;

% plot time elapse
while i<=length(Rout)
    % get Uout 
    t1 = (i-1)*options.timeStep;
    try
        minVal = inf;
        maxVal = -inf;
        for k=1:1
            Uout = interval(project(Rout{i},iDim));
            if infimum(Uout) < minVal
                minVal = infimum(Uout);
            end
            if supremum(Uout) > maxVal
                maxVal = supremum(Uout);
            end
            i = i + 1;
        end
    catch
        minVal = infimum(interval(project(Rout{i-1},iDim)));
        maxVal = supremum(interval(project(Rout{i-1},iDim)));
    end
    t2 = (i-1)*options.timeStep;
    % generate plot areas as interval hulls
    IH1 = interval([t1; minVal], [t2; maxVal]);

    plotFilled(IH1,[1 2],[0.75,0.75,0.75],'EdgeColor','none');
end

% plot simulation results
for i=1:(length(simRes.t))
    plot(simRes.t{i},simRes.x{i}*C_(iDim,:)','Color',0*[1 1 1]);
end

% plot unsafe set
unsafe = 0.00017;
plot([options.tStart,options.tFinal],[-unsafe,-unsafe],'--r');
plot([options.tStart,options.tFinal],[unsafe,unsafe],'--r');

% formatting
box on
xlabel('Time');
ylabel('y_3');
% grid on
% title('Space Station (Fixed Inputs)');

% example completed
completed = 1;




% Auxiliary functions -----------------------------------------------------

function res = specificationSafe(R)
% check if the spefication is satisfied

    res = 1;
    temp = interval(project(R,3));
    
    if supremum(temp) > 5e-4 || infimum(temp) < -5e-4 
        res = 0;
    end


%------------- END OF CODE --------------
