function example_hybrid_reach_ARCH21_heat3D_HEAT01()
% example_hybrid_reach_ARCH21_heat3D_HEAT01 -  heat conduction benchmark 
%    from the 2021 ARCH competition
%
% Syntax:  
%    example_hybrid_reach_ARCH21_heat3D_HEAT01
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%

% Author:       Niklas Kochdumper
% Written:      27-May-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


% System Dynamics ---------------------------------------------------------    

% load system matrices
load('heat5');

% construct output matric (center of block)
samples = nthroot(size(A,1),3);
indMid = ceil((size(A,1))/2);

C = zeros(1,size(A,1));
C(1,indMid) = 1;

% construct linear system object
sys = linearSys(A,B,[],C);


% Parameter ...............................................................
    
x0 = getInit(A,samples,1);
temp = diag(0.1*x0);
temp = temp(:,x0 > 0);
params.R0 = zonotope(x0,temp);

params.tFinal = 40;


% Settings ----------------------------------------------------------------

options.taylorTerms = 10;
options.zonotopeOrder = 20;
options.timeStep = 0.02;


% Reachability Analysis ---------------------------------------------------

clock = tic;
R = reach(sys,params,options);
tComp = toc(clock);


% Simulation --------------------------------------------------------------

% settings for random simulation
simOpt.points = 50;
simOpt.fracVert = 0.8;
simOpt.fracInpVert = 0.5;
simOpt.inpChanges = 20;

% simulate the system
simRes = simulateRandom(sys, params, simOpt);



% Verification ------------------------------------------------------------

goalSet = 0.10369885 + interval(-1e-4,1e-4);

% get maximum temperature
TempMax = -inf;

for i = 1:length(R.timeInterval.set)
    
    % compute maxium temperature
    val = supremum(interval(R.timeInterval.set{i}));
   
    % check if new global maximum
    if val > TempMax
        TempMax = val;
    end
end

% check if maximum temperature is in goal set
res = 0;
if in(goalSet,TempMax)
    res = 1;
end

disp(['Specification verified: ', num2str(res)]);
disp(['Computation time: ',num2str(tComp)]);


% Visualization -----------------------------------------------------------

figure; hold on; box on;

% plot reachable set
plotOverTime(R,1,'FaceColor',[.6 .6 .6],'EdgeColor','none');

% plot simulation
plotOverTime(simRes,indMid);

% formatting 
xlim([0,40]);
ylim([-0.005,0.11]);
xlabel('t');
ylabel('T');

end


% Auxiliary Functions -----------------------------------------------------

function x0 = getInit(A, samples, temp)
% returns an initial state vector

    x0 = zeros(size(A,1),1);

    if samples ~= 5 && (samples < 10 || mod(samples,10) ~= 0)
       error('Init region is not evenly divided by discretization!')
    end

    % maximum z point for initial region is 0.2 for 5 samples and 0.1 otherwise
    if samples >= 10
        max_z = (0.1/(1/samples));
    else
        max_z = (0.2/(1/samples));
    end

    for z = 0:max_z
        zoffset = z * samples * samples;

        for y = 0:(0.2/(1/samples))
            yoffset = y * samples;

            for x = 0:(0.4/(1/samples))
                index = x + yoffset + zoffset;

                x0(index+1) = temp;
            end
        end
    end
end