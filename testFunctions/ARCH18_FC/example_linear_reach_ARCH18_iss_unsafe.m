function [completed] = example_linear_reach_ARCH18_iss_unsafe(varargin)
% example_linear_reach_ARCH18_iss_unsafe - example for linear dynamics. It
%       is checked if the unsafe specification is violated as desired
%
% Syntax:  
%    example_linear_reach_ARCH18_iss_unsafe()
%
% Inputs:
%    no
%
% Outputs:
%    completed - boolean; true if the test completed without an error
%
% Example: 
%
% 
% Author:       Niklas Kochdumper
% Written:      20-June-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% load system matrices
load('iss.mat');

% construct linear system
sys = linearSys('iss',A,B,[],C);

dim=length(A);



% Options -----------------------------------------------------------------

R0 = interval(-0.0001*ones(dim,1),0.0001*ones(dim,1));
options.x0=mid(R0); %initial state for simulation

options.taylorTerms=6; %number of taylor terms for reachable sets
options.zonotopeOrder=30; %zonotope order
options.saveOrder = 1;
options.originContained=0;
options.reductionTechnique='girard';
options.linAlg = 2;
options.compOutputSet = 1;
options.outputOrder = 10;

options.verifySpecs = @specificationUnsafe;

U = interval([0;0.8;0.9],[0.1;1;1]);
options.uTrans=mid(U); %center of input set
options.U=zonotope([zeros(3,1),diag(rad(U))]); %input for reachability analysis

options.R0=zonotope(R0); %initial state for reachability analysis
options.tStart=0; %start time
options.tFinal=20; %final time
options.timeStep = 0.01;
% options.timeStep = 0.001;

% substitute options that should be optimized
if nargin == 1
    optParam = varargin{1};
    
    for i = 1:length(optParam)
        options.(optParam{i}.name) = optParam{i}.value;
    end
end




% Reachability analysis ---------------------------------------------------

% compute reachable set using zonotopes
tic
[~,~,~,Rout,~,res] = reachIss(sys,options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);

if res
   disp('Specification satisfied!'); 
else
   disp('Specification violated!'); 
end




if nargin == 0

    % Simulation --------------------------------------------------------------

    % create random simulations; RRTs would provide better results, but are
    % computationally more demanding
    runs = 100;
    fractionVertices = 1;
    fractionInputVertices = 1;
    inputChanges = 50;
    simRes = simulate_random(sys, options, runs, fractionVertices, fractionInputVertices, inputChanges);

    % simulate with RRTs
    % nrOfSamples = 100;
    % stretchingFactor = 2;
    % [~,xTraj_compl] = simulate_rrt(sys, options, nrOfSamples, Rcont, 1, stretchingFactor);

    % simulation with system eigenfrequencies as inputs
    freq = [0.7750;21.6374;9.2328];

    T = options.tStart:options.timeStep:options.tFinal;
    m = mid(U);
    r = rad(U);
    for i = 1:length(T)
        options.uTransVec(:,i) = [m(1)+r(1)*sin(freq(1)*T(i));m(2)+r(2)*sin(freq(2)*T(i));m(3)+r(3)*sin(freq(3)*T(i))];
    end

    runs = 1;
    fractionVertices = 1;
    fractionInputVertices = 1;
    inputChanges = 1;
    simEig = simulate_random(sys, options, runs, fractionVertices, fractionInputVertices, inputChanges);




    % Visualization -----------------------------------------------------------

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
            for k=1:3
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

        plotFilled(IH1,[1 2],[.75 .75 .75],'EdgeColor','none');
    end


    % plot simulation results RRT
    % t = options.tStart;
    % for i = 1:length(xTraj_compl)
    %    temp = xTraj_compl{i};
    %    for j = 1:length(temp)
    %       len = size(temp{j},1);
    %       T = linspace(t,t+options.timeStep,len);
    %       plot(T,C(iDim,:)*temp{j}','b');
    %    end
    %    t = t + options.timeStep;
    % end

    % plot simulation results
    for i=1:(length(simRes.t))
        plot(simRes.t{i},simRes.x{i}*C(iDim,:)','Color',0*[1 1 1]);
    end

    % plot simulation results (eigenfrequencies)
    for i=1:(length(simEig.t))
        plot(simEig.t{i},simEig.x{i}*C(iDim,:)','k');
    end

    % plot unsafe set
    unsafe = 0.0005;
    plot([options.tStart,options.tFinal],[-unsafe,-unsafe],'--r');
    plot([options.tStart,options.tFinal],[unsafe,unsafe],'--r');

    % formatting
    box on
    xlabel('Time');
    ylabel('y_3');
    grid on
    % title('Space Station (Time-varying Inputs)');
    %grid on

end



% Output arguments --------------------------------------------------------

completed = 1;



% Auxiliary functions -----------------------------------------------------

function res = specificationUnsafe(R)
% check if the spefication is satisfied

    res = 1;
    temp = interval(project(R,3));
    
    if supremum(temp) > 5e-4 || infimum(temp) < -5e-4 
        res = 0;
    end


%------------- END OF CODE --------------
