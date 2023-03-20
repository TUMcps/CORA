function example_linear_reach_ARCH22_heat3D_HEAT02()
% example_linear_reach_ARCH22_heat3D_HEAT02 - heat conduction benchmark 
%    from the 2022 ARCH competition
%
% Syntax:  
%    example_linear_reach_ARCH22_heat3D_HEAT02
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      02-June-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


% System Dynamics ---------------------------------------------------------    

% load system matrices
load('heat10');

% construct output matric (center of block)
samples = nthroot(size(A,1),3);
indMid = 556;

C = zeros(1,size(A,1));
C(1,indMid) = 1;


% Parameters --------------------------------------------------------------
    
x0 = getInit(A,samples,1);
temp = diag(0.1*x0);
temp = temp(:,x0 > 0);
R0sim = zonotope(x0,temp);

params.tFinal = 40;

% construct transpose linear system object
sys = linearSys(A,0,[],B);
sysTrans = linearSys(A.',0,[],R0sim.Z.');

% correct input
params.R0 = zonotope(C.');


% Settings ----------------------------------------------------------------

options.taylorTerms = 10;
options.zonotopeOrder = 2;
options.timeStep = 0.05;
options.linAlg = 'krylov';
options.krylovError = eps;
options.krylovStep = 20;


% Reachability Analysis ---------------------------------------------------

clock = tic;
Rcont = reach(sysTrans,params,options);
tComp = toc(clock);


% Simulation --------------------------------------------------------------

% settings for random simulation
simOpt.points = 1;
simOpt.fracInpVert = 0.5;

params.R0 = R0sim;

% simulate the system
simRes = simulateRandom(sys, params, simOpt);



% Verification ------------------------------------------------------------

goalSet = 0.02966356 + interval(-1e-4,1e-4);

% get maximum temperature
TempMax = -inf;

for i = 1:length(Rcont(1).timeInterval.set)
    
    % compute maxium temperature
    Rtmp = Rcont(1).timeInterval.set{i}.Z;
    Rtmp = [sum(Rtmp(1,:)); reshape(Rtmp(2:end,:),[],1)];
    Rset{i} = zonotope(Rtmp.');
    val = supremum(interval(Rset{i}));
   
    % check if new global maximum
    if val(1) > TempMax
        TempMax = val;
    end
end

% check if maximum temperature is in goal set
res = 0;
if contains(goalSet,TempMax)
    res = 1;
end

disp(['Specification verified: ', num2str(res)]);
disp(['Computation time: ',num2str(tComp)]);


% Visualization -----------------------------------------------------------

figure; hold on; box on;

% plot reachable set
t = 0;

for i = 1:length(Rset)
    
    % get intervals
    intX = interval(Rset{i});
    intT = interval(t,t+options.timeStep);
    
    int = cartProd(intT,intX(1));
    
    % plot the interval
    plot(int,[1,2],'FaceColor',CORAcolor("CORA:reachSet"));
    
    % update time
    t = t + options.timeStep;
end

% plot simulation
for i = 1:length(simRes.x)
   plot(simRes.t{i},simRes.x{i}(:,indMid),'Color',CORAcolor("CORA:simulations"));
end

% formatting 
xlim([0,40]);
ylim([-0.005,0.035]);
xlabel('t');
ylabel('T');


end


% Auxiliary Functions -----------------------------------------------------

function x0 = getInit(A, samples, temp)
% returns an initial state vector

    x0 = zeros(size(A,1),1);

    if samples ~= 5 && (samples < 10 || mod(samples,10) ~= 0)
        throw(CORAerror('CORA:specialError',...
            'Init region is not evenly divided by discretization!'));
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

%------------- END OF CODE --------------