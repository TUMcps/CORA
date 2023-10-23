function completed = example_hybrid_reach_03_PLLnoSat()
% example_hybrid_reach_03_PLLnoSat - example for hybrid dynamics; this 
%    example uses continuization to speed up the verification process;
%    this example can be found in [1] and [2].
%
% Syntax:
%    completed = example_hybrid_reach_03_PLLnoSat()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false
%
% References:
%    [1] Althoff, M.; Rajhans, A.; Krogh, B. H.; Yaldiz, S.; Li, X. & Pileggi, L. 
%        Formal Verification of Phase-Locked Loops Using Reachability
%        Analysis and Continuization. Proc. of the Int. Conference on
%        Computer Aided Design, 2011, 659-666.
%    [2] Althoff, M.; Rajhans, A.; Krogh, B. H.; Yaldiz, S.; Li, X. & Pileggi, L. 
%        Formal Verification of Phase-Locked Loops Using Reachability Analysis 
%        and Continuization Communications of the ACM, 2013, 56, 97-104

% Authors:       Matthias Althoff
% Written:       20-January-2010
% Last update:   26-May-2018
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%initial state
minPhase = -0.5;
maxPhase = -0.4;
R0ih=interval([0.34; -0.01; -0.01; minPhase], [0.36; 0.01; 0.01; maxPhase]);

% instantiate hybrid automaton for simulation
[HAsim, A, U, U_delay, Uloc, c] = PLL();


% Parameters --------------------------------------------------------------

params.R0 = zonotope(R0ih);
params.U = U;
params.U_delay = U_delay;
params.Uloc = Uloc;
params.startLoc = 2; %initial location 
params.uLoc{1} = center(params.Uloc{1});
params.uLoc{2} = randPoint(params.Uloc{2});
params.uLoc{3} = randPoint(params.Uloc{3});
uTmp = center(params.Uloc{1});
uTmp(1:3) = params.uLoc{2}(1:3) + params.uLoc{3}(1:3);
params.uLoc{4} = uTmp;
params.tStart = 0;
params.tFinal = Inf;


% Reachability Settings ---------------------------------------------------

options.taylorTerms = 10;
%options.zonotopeOrder = 5; 
options.Rlock=interval([-100; -100; -100; -1/3600],[100; 100; 100; 1/3600]);               


% Reachability Analysis ---------------------------------------------------

% no saturation
tic
R = aux_reachPLL_general_noSat(A, c, params, options);
tComp_reach = toc;


% Simulation --------------------------------------------------------------

plotRange = 1:200;
%simCycles = length(R);
simCycles = plotRange(end);

%obtain random simulation results
%samples = 30;
samples = 3;

tic
for iSim = 1:samples
    % set initial state, input
    if iSim<=15
        p = randPoint(params.R0,1,'extreme');
        params.x0 = [p;0;0];                % initial state for simulation
    else
        p = randPoint(params.R0);
        params.x0 = [p;0;0];                % initial state for simulation
    end
    
    if params.x0(4)~=params.x0(5)
        
        %set initial location
        if params.x0(4)<params.x0(5)
            params.loc = 2; %initial location 
        else
            params.loc = 1; %initial location 
        end
    
        %simulate hybrid automaton
        x_cycle = aux_simulatePLL(HAsim,params,simCycles); 

        %get simulation results
        simRes{iSim} = x_cycle;
    end
    
    iSim
end

tComp_sim = toc;


% Visualization -----------------------------------------------------------

% projections for plotting
Psim = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 -1 0];
dims = [1 2; 2 3; 1 4; 3 4];

% loop over all projections
for i = 1:4
    
    figure; hold on;
    
    %plot reachable set
    for iSet = plotRange
        Rproj = reduce(project(R{iSet},dims(i,:)),'girard',20);
        plot(Rproj,[1 2],'FaceColor',colorblind('b'));
    end
    
    %plot simulation results
    for iSet = plotRange
        for n = plotRange
            %plot results
            try
                xProj = Psim*simRes{n}(iSet,:)';
                plot(xProj(dims(i,1)),xProj(dims(i,2)),...
                    'Color',colorblind('y'),'Marker','.');
                %plot(simRes{n}(iSet,dims(i,1)),simRes{n}(iSet,dims(i,2)),'ro');
            catch
            end
        end
    end
end

%example completed
completed = true;

end


% Auxiliary functions -----------------------------------------------------

% main function for reachability analysis
function Rout = aux_reachPLL_general_noSat(A_tmp, c_tmp, params, options)

% projection to first four dimensions
A = A_tmp(1:4, 1:4);
c = c_tmp(1:4);
U = project(params.U, 1:4);
U_delay = project(params.U_delay, 1:4);

%obtain variables
R{1} = params.R0; 

%compute phase and voltage velocity boundaries
viInt = interval(0,0.7);
vpInt = interval(-4,12);
[vMinPhase,vMaxPhase]=aux_phaseVelBound(A, c, viInt, vpInt);

%set cycle and delay time
t_cycle = 1/27;
t_delay = 50e-6;

%compute interval hull
iCycle = 1;

%precompute constant matrices
Apower = aux_powers(A,30);
prev_taylorterms = options.taylorTerms;
taylorTerms = 30;
eAtInt_delay = aux_inputExponential(A, Apower, t_delay, taylorTerms);
Rconst = aux_inputExponential(A, Apower,t_cycle,taylorTerms)*zonotope(c);
taylorTerms = prev_taylorterms;

notLocked = 1;
lockingViolation = 1;

while notLocked

    %compute time limits when charge pump is on
    [t_min,t_max] = aux_timeBound_phase(R{iCycle},vMinPhase,vMaxPhase);
    
    %determine delta T
    deltaT = t_max - t_min;
    
    %compute phase velocity interval
    vInt = interval(vMinPhase,vMaxPhase);
    
    %compute time delay solution
    eAt_unc = aux_uncertainTimeExponential(A, Apower, deltaT, taylorTerms);
    Rdelay = expm(A*(t_cycle-t_max-t_delay)) * eAt_unc * eAtInt_delay* U_delay;
    
    %compute input matrix for t=[t_min, t_max] 
    matZinput = aux_inputMatrix(Apower, deltaT, U, vInt, taylorTerms);
    
    if t_min>0
    
        %compute input solution for t=[0,t_min] 
        eAtInt_t_min = aux_inputExponential(A, Apower, t_min, taylorTerms);
        Rinput = eAtInt_t_min * (U + c);

        %complete solution
        prev_taylorterms = taylorTerms;
        taylorTerms = 30;
        Rconst_2 = aux_inputExponential(A, Apower,(t_cycle-t_min),taylorTerms)*zonotope(c);
        taylorTerms = prev_taylorterms;

        R_t_min = expm(A*t_min)*R{iCycle} + Rinput;
        R_new = (expm(A*(t_cycle-t_min))+matZinput)*R_t_min + Rdelay + Rconst_2; 
        %Rout
        Rout{iCycle} = R{iCycle};
    else
        %complete solution
        R_new = (expm(A*t_cycle)+expm(A*(t_cycle-t_max))*matZinput)*R{iCycle} + Rdelay + Rconst; 
        %Rout
        Rout{iCycle} = R{iCycle} + matZinput*R{iCycle};
    end
    
    %reset 
    R_new = R_new + [0;0;0;-1];
    
    %reduce 
    R{iCycle+1} = reduce(R_new,'girard',100);  
    
    %interval of previous reachable set
    IH=interval(Rout{iCycle});
    if iCycle>200
        IHold = interval(Rout{iCycle-200});
    else
        IHold = IH;
    end
    
    %compute if in locking
    inLocking = (IH<=options.Rlock);
    inLockingOld = (IHold<=options.Rlock);
    
    if inLocking && inLockingOld
        
        %simplify to box if there was a locking violation
        if lockingViolation
            %box enclosure
            R{iCycle+1} = zonotope(IH);
            IHstart = IH;
            %update locking violation
            lockingViolation = 0;
        end
        
        %check enclosure
        if (IH<IHstart) && ~lockingViolation
            notLocked = 0;
        end
    else
        lockingViolation =1;
    end
    
    iCycle = iCycle+1;
end

end


%compute phase velocity boundaries
function [vMin,vMax]=aux_phaseVelBound(A,c,vpInt,viInt)

vInt = A(4,:)*[viInt; 0; vpInt; 0] + c(4);

vMin = infimum(vInt);
vMax = supremum(vInt);

end


%compute time intervals for charge pump on and off
function [t_min,t_max]=aux_timeBound_phase(R,vMin,vMax)

%obtain range of Phi_v
Phi_IH = interval([0 0 0 1]*R);
PhiMin = min(abs(infimum(Phi_IH)),abs(supremum(Phi_IH)));
PhiMax = max(abs(infimum(Phi_IH)),abs(supremum(Phi_IH)));

%t_on
if infimum(Phi_IH)*supremum(Phi_IH)>0
    t_min = PhiMin/vMax;
else
    t_min = 0;
end
t_max = PhiMax/vMin;

end


function eAtInt = aux_inputExponential(A,Apower,r,taylorTerms)

%compute E
E = aux_remainder(A,r,taylorTerms);

dim = length(Apower{1});
Asum = r*eye(dim);
%compute higher order terms
for i=1:taylorTerms
    %compute factor
    factor = r^(i+1)/factorial(i+1);    
    %compute sums
    Asum = Asum + Apower{i}*factor;
end

%compute exponential due to constant input
eAtInt = Asum + E*r;

end


function eAt = aux_uncertainTimeExponential(A,Apower,deltaT,taylorTerms)

%compute Apowers
E = aux_remainder(A,deltaT,taylorTerms);

%define time interval as matrix zonotope
for i=1:taylorTerms
    gen_tmp{1} = 0.5*deltaT^i;
    t_int{i} = matZonotope(0.5*deltaT^i,gen_tmp); 
end

dim = length(Apower{1});
gen{1} = zeros(dim);
Asum = matZonotope(eye(dim),gen);
%compute higher order terms
for i=1:taylorTerms
    %compute factor
    factor = t_int{i}*(1/factorial(i));    
    %compute sums
    Asum = Asum + Apower{i}*factor;
end

%compute exponential due to constant input
eAt = intervalMatrix(Asum) + E;

end


function Apower = aux_powers(A,taylorTerms)

%initialize 
Apower = cell(1,taylorTerms+1);
Apower{1} = A;  
    
%compute powers for each term and sum of these
for i=1:taylorTerms
    %compute powers
    Apower{i+1}=Apower{i}*A;
end   

end

function E = aux_remainder(A,r,taylorTerms)

%compute absolute value bound
M = abs(A*r);
dim = length(M);

%compute exponential matrix
eM = expm(M);

%compute first Taylor terms
Mpow = eye(dim);
eMpartial = eye(dim);
for i=1:taylorTerms
    Mpow = M*Mpow;
    eMpartial = eMpartial + Mpow/factorial(i);
end

W = eM-eMpartial;

%instantiate remainder
E = intervalMatrix(zeros(dim),W);

end

function E = aux_expMatRemainder(Apower, t, taylorTerms)

%compute auxiliary matrix
tmpMat = 1/factorial(4)*(t^4)*Apower{4};

%compute E_1 center, generators
E_center = 0.5*tmpMat;
E_gen{1} = 0.5*tmpMat;

for i=5:taylorTerms
    tmpMat = 0.5*(t^i)*Apower{i};
    
    %compute E center, generators
    E_center = E_center + 1/factorial(i)*tmpMat;
    E_gen{end+1} = 1/factorial(i)*tmpMat;
end

%instantiate matrix zonotopes
E_zono = matZonotope(E_center, E_gen);

%compute remainders
Erem = exponentialRemainder(intervalMatrix(0.5*Apower{1}*t,0.5*Apower{1}*t),taylorTerms);

%final result
E = intervalMatrix(E_zono) + Erem;

end

function matZinput = aux_inputMatrix(Apower, Tmax, U, vInt, taylorTerms)

%initialize
dim = length(Apower{1});

t_intMat = intervalMatrix(0.5*Tmax,0.5*Tmax); %other time delay: 50e-3

T(1) = interval(1/2*Tmax,Tmax); %other time delay: 50e-3
T(2) = interval(1/6*Tmax^2,1/2*Tmax^2); %other time delay: 50e-3
T(3) = interval(1/8*Tmax^3,1/6*Tmax^3); %other time delay: 50e-3

%convert to matrix zonotopes
for i=1:length(T)
    gen{1} = rad(T(i));
    t_zonMat{i} = matZonotope(center(T(i)),gen);
end

%compute matrix zonotopes
E = aux_expMatRemainder(Apower, Tmax, taylorTerms);

%compute factor
factor = (1/vInt);
matFactor_int = intervalMatrix(center(factor),rad(factor));
matFactor_zono = matZonotope(matFactor_int);

%zonotope part of input matrix
inputMat_zono = matFactor_zono*(eye(dim) + t_zonMat{1}*Apower{1}...
                + t_zonMat{2}*Apower{2} + t_zonMat{3}*Apower{3});        
inputMat_int = matFactor_int*E*(eye(dim) + 1/factorial(2)*t_intMat*Apower{1} ...
                + 1/factorial(3)*t_intMat^2*Apower{2} + 1/factorial(4)*t_intMat^3*Apower{3} + E);
%----------------------------------------------------------------------

%INPUT-----------------------------------------------------------------


%compute fourth column
U = interval(U);
U_inf = infimum(U);
U_sup = supremum(U);
U_intMat = intervalMatrix(0.5*(U_inf+U_sup), 0.5*(U_sup-U_inf));
fourthCol = (-1)*(intervalMatrix(inputMat_zono) + inputMat_int)*U_intMat;
fourthCol_center = center(fourthCol.int);
fourthCol_gen = rad(fourthCol.int);

Acenter = zeros(dim);
Adelta = zeros(dim);
Acenter(:,4) = fourthCol_center;
Adelta(:,4) = fourthCol_gen;
matZinput = intervalMatrix(Acenter,Adelta);

end

function x_cycle = aux_simulatePLL(obj,params,simCycles)

%initialize intermediate state
x_cycle=params.x0'; %intermediate state at transitions

%while final time not reached
while length(x_cycle(:,1)) < simCycles

    %simulate within the actual location
    params.u = params.uLoc{params.loc};
    [~,~,params.loc,params.x0]=simulate(obj.location(params.loc),params);
    
    %store results at the end of location 1
    if params.loc==2 && (params.x0(4)<=params.x0(5))
        x_cycle(end+1,:) = params.x0';
    end
    if params.loc==3 && (params.x0(4)>=params.x0(5)) % loc 4 was used previously; check
        x_cycle(end+1,:) = params.x0';
    end
end

end

% ------------------------------ END OF CODE ------------------------------
