function [Rout,IHstart] = reachPLL_general_sat(obj, nonLinDyn, options)
% reachPLL_general_sat - computes the reachable set of the PLL considering
% saturation
%
% Syntax:  
%    R = reachPLL_general_sat(obj, options)
%
% Inputs:
%    options - options for reachability analysis of the system
%
% Outputs:
%    R - reachable set
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      14-November-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


%obtain variables
A = options.sys.A{1};
U = options.U;
U_delay = options.U_delay;
U_full = options.U_full;
c = options.sys.c{1};
Rcycle{1} = options.R0; 
uTmp = interval(U);
u_max = uTmp(:,2);
u_mid = center(U);

viInt = options.viInt;
vpInt = options.vpInt;


%compute phase and voltage velocity boundaries
[vMinPhase,vMaxPhase] = phaseVelBound(A, c, viInt, vpInt);

%set cycle and delay time
t_cycle = 1/27;
t_delay = 50e-6;

%compute interval hull
iCycle = 1;
iLock = 0;


%precompute constant matrices
Apower = powers(A,30);
prev_taylorterms = options.taylorTerms;
options.taylorTerms = 30;
eAtInt_delay = inputExponential(A, Apower, t_delay, options);
Rconst = inputExponential(A, Apower,t_cycle,options)*zonotope(c);
options.taylorTerms = prev_taylorterms;

negPhaseDiff = 1;
notLocked = 1;
lockingViolation = 1;
iCycle = 1;

while negPhaseDiff
    %PHASE 1---------------------------------------------------------------
    %compute time limits when charge pump is on
    [t_min,t_max] = timeBound_phase(Rcycle{iCycle},vMinPhase,vMaxPhase);
    
    %determine delta T
    deltaT = t_max - t_min;
    
    %compute phase velocity interval
    vInt = infsup(vMinPhase,vMaxPhase);
    
    %reachable set until t_min
    %set time step
    options.timeStep = t_min/ceil(t_min/options.timeStep);
    options.U = U;
    options.R0 = Rcycle{iCycle};
    Rlast = nonlinReach(nonLinDyn, options, t_min); %return time point not time interval reachable set!
    
    %reachable set from t_min to t_max
    options.timeStep = deltaT/ceil(deltaT/options.timeStep);
    options.U = U_full;
    options.R0 = Rlast;
    R = nonlinReach(nonLinDyn, options, t_max);
    
    %reachable set from t_max to t_max+t_delay
    options.timeStep = t_delay;
    options.U = U_delay;
    options.R0 = Rlast;
    R = nonlinReach(nonLinDyn, options, t_max + t_delay);
    
    %compute solution from t_max to t_cycle
    iCycle = iCycle + 1;
    Rcycle{iCycle} = expm(A*(t_cycle - t_max - t_delay))*R{end}{1};
end

while notLocked

    %PHASE 1---------------------------------------------------------------
    %compute time limits when charge pump is on
    [t_min,t_max] = timeBound_phase(R{iCycle},vMinPhase,vMaxPhase);
    
    %determine delta T
    deltaT = t_max - t_min;
    
    %compute phase velocity interval
    vInt = infsup(vMinPhase,vMaxPhase);
    
    %compute worst case CP values
    RfullInput = inputExponential(A, Apower, t_max, options)*zonotope(u_max);
    RmaxCP = expm(A*t_max)*R{iCycle};
    IHmaxCP = interval(RmaxCP + RfullInput);
    
    %linearize
    [A, B, Ulin] = linearizePll(obj.mFile, U, t_max);
    
    %compute time delay solution
    eAt_unc = uncertainTimeExponential(A, Apower, deltaT, options);
    Rdelay = expm(A*(t_cycle-t_max-t_delay)) * eAt_unc * eAtInt_delay* U_delay;
    
    %compute input matrix for t=[t_min, t_max] 
    matZinput = inputMatrix(Apower, deltaT, U, vInt, options);
    
    if t_min>0
    
        %compute input solution for t=[0,t_min] 
        eAtInt_t_min = inputExponential(A, Apower, t_min, options);
        Rinput = eAtInt_t_min * (U + c);

        %complete solution
        prev_taylorterms = options.taylorTerms;
        options.taylorTerms = 30;
        Rconst_2 = inputExponential(A, Apower,(t_cycle-t_min),options)*zonotope(c);
        options.taylorTerms = prev_taylorterms;

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
    
    %FOR COMPARISON:
%     %compute input solution for t=[t_min,t_max] 
%     eAtInt = inputExponential(A, Apower,deltaT,options);
%     Rinput_test = expm(A*(t_cycle-t_max)) * eAtInt * U;
%     Rinput_test2= matZinput*R_t_min;
    
    %reduce 
    R{iCycle+1} = reduce(R_new,'girard',100);
   
    
    %interval of previous reachable set
    IH=interval(Rout{iCycle});
    if iCycle>200
        IHold = interval(Rout{iCycle-200});
    else
        IHold = IH;
    end
   
    if mod(iCycle,500)==0
        disp('next 500 cycles')
        IH
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
            
            Rold=Rout{iCycle-200};
            Rstart=Rout{iCycle};
        end
        
        P=[0 0 0 1; 1 0 0 0];
        Rred = reduce(P*Rout{iCycle},'girard',10);
        plot(Rred,[1 2],'lightgray');
        
        %check enclosure
        if (IH<IHstart) && ~lockingViolation
            plot(options.Rlock,[4 1]);
            plot(Rold,[4 1]);
            plot(Rstart,[4 1]);
            plot(IHstart,[4 1]);
            plot(Rout{iCycle},[4 1]);
            
            axis([-0.0005 0.0005 0.349 0.351]);
            notLocked = 0;
        end
    else
        lockingViolation =1;
    end
    
    iCycle = iCycle+1;
end


 
%compute phase velocity boundaries
function [vMin,vMax]=phaseVelBound(A,c,vpInt,viInt)

vInt = A(4,:)*[viInt; 0; vpInt; 0] + c(4);

vMin = inf(vInt);
vMax = sup(vInt);



%compute time intervals for charge pump on and off
function [t_min,t_max]=timeBound_phase(R,vMin,vMax)

%obtain range of Phi_v
Phi_IH = interval([0 0 0 1]*R);
PhiMin = min(abs(Phi_IH(:,1)),abs(Phi_IH(:,2)));
PhiMax = max(abs(Phi_IH(:,1)),abs(Phi_IH(:,2)));

%t_on
if Phi_IH(:,1)*Phi_IH(:,2)>0
    t_min = PhiMin/vMax;
else
    t_min = 0;
end
t_max = PhiMax/vMin;


function eAtInt = inputExponential(A,Apower,r,options)

%compute E
E = remainder(A,r,options);

dim = length(Apower{1});
Asum = r*eye(dim);
%compute higher order terms
for i=1:options.taylorTerms
    %compute factor
    factor = r^(i+1)/factorial(i+1);    
    %compute sums
    Asum = Asum + Apower{i}*factor;
end

%compute exponential due to constant input
eAtInt = Asum + E*r;


function eAt = uncertainTimeExponential(A,Apower,deltaT,options)

%compute Apowers
E = remainder(A,deltaT,options);

%define time interval as matrix zonotope
for i=1:options.taylorTerms
    gen_tmp{1} = 0.5*deltaT^i;
    t_int{i} = matZonotope(0.5*deltaT^i,gen_tmp); 
end

dim = length(Apower{1});
gen{1} = zeros(dim);
Asum = matZonotope(eye(dim),gen);
%compute higher order terms
for i=1:options.taylorTerms
    %compute factor
    factor = t_int{i}*(1/factorial(i));    
    %compute sums
    Asum = Asum + Apower{i}*factor;
end

%compute exponential due to constant input
eAt = intervalMatrix(Asum) + E;



function Apower = powers(A,taylorTerms)

%initialize 
Apower = cell(1,taylorTerms+1);
Apower{1} = A;  
    
%compute powers for each term and sum of these
for i=1:taylorTerms
    %compute powers
    Apower{i+1}=Apower{i}*A;
end   


function E = remainder(A,r,options)

%compute absolute value bound
M = abs(A*r);
dim = length(M);

%compute exponential matrix
eM = expm(M);

%compute first Taylor terms
Mpow = eye(dim);
eMpartial = eye(dim);
for i=1:options.taylorTerms
    Mpow = M*Mpow;
    eMpartial = eMpartial + Mpow/factorial(i);
end

W = eM-eMpartial;

%instantiate remainder
E = intervalMatrix(zeros(dim),W);


function [E] = expMatRemainder(Apower, t, options)

%compute auxiliary matrix
tmpMat = 1/factorial(4)*(t^4)*Apower{4};

%compute E_1 center, generators
E_center = 0.5*tmpMat;
E_gen{1} = 0.5*tmpMat;

for i=5:options.taylorTerms
    tmpMat = 0.5*(t^i)*Apower{i};
    
    %compute E center, generators
    E_center = E_center + 1/factorial(i)*tmpMat;
    E_gen{end+1} = 1/factorial(i)*tmpMat;
end

%instantiate matrix zonotopes
E_zono = matZonotope(E_center, E_gen);

%compute remainders
Erem = exponentialRemainder(intervalMatrix(0.5*Apower{1}*t,0.5*Apower{1}*t),options.taylorTerms);

%final result
E = intervalMatrix(E_zono) + Erem;


function matZinput = inputMatrix(Apower, Tmax, U, vInt, options)

%initialize
dim = length(Apower{1});

t_intMat = intervalMatrix(0.5*Tmax,0.5*Tmax); %other time delay: 50e-3

T(1) = infsup(1/2*Tmax,Tmax); %other time delay: 50e-3
T(2) = infsup(1/6*Tmax^2,1/2*Tmax^2); %other time delay: 50e-3
T(3) = infsup(1/8*Tmax^3,1/6*Tmax^3); %other time delay: 50e-3

%convert to matrix zonotopes
for i=1:length(T)
    gen{1} = rad(T(i));
    t_zonMat{i} = matZonotope(mid(T(i)),gen);
end

%compute matrix zonotopes
E = expMatRemainder(Apower, Tmax, options);

%compute factor
factor = (1/vInt);
matFactor_int = intervalMatrix(mid(factor),rad(factor));
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
U_inf = U(:,1);
U_sup = U(:,2);
U_intMat = intervalMatrix(0.5*(U_inf+U_sup), 0.5*(U_sup-U_inf));
fourthCol = (-1)*(intervalMatrix(inputMat_zono) + inputMat_int)*U_intMat;
fourthCol_center = mid(fourthCol.int);
fourthCol_gen = rad(fourthCol.int);

Acenter = zeros(dim);
Adelta = zeros(dim);
Acenter(:,4) = fourthCol_center;
Adelta(:,4) = fourthCol_gen;
matZinput = intervalMatrix(Acenter,Adelta);


function [Alin, Blin, xlin] = linearizePll(contDyn, Rinit, u_mid, t_max, options)

%simulate
%generate overall options
stepsizeOptions = odeset('MaxStep',0.2*t_max);
opt = odeset(stepsizeOptions);

%simulate for half the time
options.u = u_mid;
x0 = center(Rinit);
[contDyn,t,x] = simulate(contDyn,options,0,0.5*t_max,x0,opt); %linearize at 0.5t_max

%substitute p into the Jacobian with respect to x and u to obtain the
%system matrix A and the input matrix B
xlin = x(end,:)';
[Alin, Blin] = jacobian(xlin,u_mid);


%compute worst case CP values
function [vi_min, vp1_min, vi_max, vp1_max] = worstCaseCPvalues(contDyn, Rinit, u_max, t_max, options)

%simulate
%generate overall options
stepsizeOptions = odeset('MaxStep',0.2*t_max);
opt = odeset(stepsizeOptions);

%simulate 
options.u = u_max;
IH = interval(Rinit);
x_min = IH(:,1);
x_max = IH(:,2);

if x_max(4)<0
    [contDyn,t,x] = simulate(contDyn,options,0,t_max,x_max,opt); 

    %extract values
    vi_min = x_min(1);
    vp1_min = x_min(2);
    vi_max = x(end,1);
    vp1_max = x(end,2);
end

function Rlast = nonlinReach(sys, options, tFinal)

%initialize reachable set computations
[sys, Rnext, options] = initReach(sys, options.R0, options);

%while final time is not reached
t=options.tStart;
iSet=1;

while t<tFinal
    
    %increment time and set counter
    t = t+options.timeStep;
    
    %compute next reachable set
    [Rnext,options]=reach(sys,Rnext,options);
end

Rlast = Rnext.tp{1};

%------------- END OF CODE --------------