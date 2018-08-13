function [R, IHintersectAll, timeStep_max, extraCycles] = reachPLL_general(obj, timeStep_max_prev, options)
% reachPLL_general - computes the reachable set of the PLL 
%
% Syntax:  
%    [R] = reachPLL_lock(obj, options)
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
% Written:      14-December-2010
% Last update:  22-December-2010
%               12-January-2010
% Last revision:---

%------------- BEGIN CODE --------------

%define A and A_input!!

%obtain variables
initSatInd = options.initSatInd;

A = options.sys.A{1};
A_input = options.sys.A{initSatInd};
c = options.sys.c{initSatInd};
U = options.U{initSatInd};
guard = options.Rguard{initSatInd};
inv = options.Rinv{initSatInd};
R{1} = options.R0; 

dim = length(A);

IHintersectAll = interval();
timeStep_max=0;
extraCycles=0;

%generate Uuncertain
U_ih = interval(U) + (-c);
Uuncertain = zonotope(interval(0*U_ih(:,2), abs(U_ih(:,2)))) + c;

%instantiate linear system
linSys = linearSys('linearInputDynamics',A_input,eye(dim)); 

%options for reachability analysis
options.uTrans = center(U);
options.U = U + (-options.uTrans);

if options.inputPhaseInit
    %compute initial set
    [linSys,Rfirst,options] = initReach(linSys,R{1},options);

    %init
    IH = interval();
    t = 0;
    Rnext = Rfirst.ti;
    IH = interval(Rfirst.ti);
    while t<timeStep_max_prev
        %compute next reachable set
        [Rnext, options] = reach(linSys, Rnext, options);

        %boxing
        IHnew = interval(Rnext.ti);

        %unify boxes
        IH = IH | IHnew;

        %increment time 
        t = t+options.timeStep;
    end
    %overwrite initial set
    R{1} = zonotope(IH);
    R{1} = zonotopeBundle(R);
    %intersect with invariant
    R{1} = R{1} & zonotope(inv);
    R{1} = shrinkIH(R{1});
end


%compute phase and voltage velocity boundaries
%set viInt and vpInt
viInt = infsup(-0.5,0.5);
vpInt = infsup(-0.5,0.5);
vSat = 0.35; %<-- change vSat here!
[vMinPhase,vMaxPhase]=phaseVelBound(A_input, c, viInt, vpInt);
[vMinVoltage,vMaxVoltage] = voltageVelBound(A_input,vpInt,U);
D = derivativeSet(A, U);

%compute cycle time
t_cycle = 1/27;

%compute interval hull
IH = interval(R{1});


% %init
% Vall = cell(1,length(options.Pinv));
% for iSet = 1:length(options.Pinv);
%     Vall{iSet} = [];
% end

iCycle=1;

while ~isempty(IH) && ~(IH<=options.Rgoal{1}) && ~(IH<=options.Rgoal{2})

    %PHASE 1---------------------------------------------------------------
    %compute time limits when charge pump is on
    [t_min,t_max] = timeBound_phase(R{iCycle},vMinPhase,vMaxPhase);
    t_on_voltage = timeBound_voltage_conservative(R{iCycle},vMaxVoltage,vSat);
    
    %determine delta T
    deltaT = t_max - t_min;
    
%     %determine saturation
%     if (t_on_voltage(1)>=t_max) && (t_on_voltage(2)>=t_max)
%         %no saturation
%         satInd = 1;
%     elseif (t_on_voltage(1)<t_max) && (t_on_voltage(2)>=t_max)
%         %integral saturation
%         satInd = 2;
%     elseif (t_on_voltage(1)>=t_max) && (t_on_voltage(2)<t_max)
%         %proportional saturation
%         satInd = 3;
%     else
%         %integral and proportional saturation
%         satInd = 4;
%     end
    
    %copute reachable set of a single cycle
    R_new = singleCycleReach(A, A_input, U, Uuncertain, t_min, t_max, t_cycle, R{iCycle}, options);
   
    
    %compute solution for intersection with voltage bounds
    satReached = 0;
    if initSatInd==1 && t_on_voltage(2) < t_max
        satReached = 1;
    elseif initSatInd==3 && t_on_voltage(1) < t_max
        satReached = 1;
    end
    
    if satReached
        
        if t_max > timeStep_max
            timeStep_max = t_max;
        end

        %compute initial set
        [linSys,Rfirst,options] = initReach(linSys,R{iCycle},options);
        
        %init
        inInv = 1;
        IHintersect = interval();
        t = 0;
        Rnext = Rfirst.ti;
        while inInv && (t<t_max)
            %compute next reachable set
            [Rnext, options] = reach(linSys, Rnext, options);
            
            %intersection
            IHintersectNew = interval(Rnext.ti & guard);
            
            %unify intersections
            if ~isempty(IHintersectNew)
                if isempty(IHintersect)
                    IHintersect = IHintersectNew;
                else
                    IHintersect = IHintersect | IHintersectNew;
                end
            end
            
            %increment time and set counter
            t = t+options.timeStep;

            inInv=in(inv,Rnext.ti); %check if reachable set is in invariant 
        end
 

        %store result
        if isempty(IHintersectAll)
            IHintersectAll = IHintersect;
        else
            IHintersectAll = IHintersectAll | IHintersect;
        end
    end
    
    
    %intersect with invariant
    R_new = R_new & zonotope(inv);
    R_new = shrinkIH(R_new);
    %IH_intersect = interval(R_intersect);
    
    if ~isempty(R_new)
        IH=interval(R_new);
        R{iCycle+1} = reduce(R_new,'girard',60);
    else
        %initialize empty interval
        IH = interval();
    end
    
    %INTRODUCE PURE SPLITTING HERE...
    
    %if mod(iCycle,400)==0
    if mod(iCycle,10)==0
        %obtain linear map
        t = iCycle*t_cycle;
        linMap = expm(A*t) + 0.01*eye(dim);
        %simplify
        R{iCycle+1} = simplifySet(R{iCycle+1}, linMap, A, 1);
        
        %contract if not in saturation
        if satReached==0
            %one cycle phase difference approximation
            centerDiff_cycle = center(R{iCycle+1}.Z{2}) - center(R{iCycle}.Z{2});
            phaseDiff_cycle = centerDiff_cycle(4);
            
            %split reachable set in x_4 direction
            IHmat = get(IH,'intervals');
            
            IHmat1 = IHmat;
            IHmat1(4,2) = 0.5*(IHmat(4,2) + IHmat(4,1));
            IH1 = interval(IHmat1(:,1), IHmat1(:,2));
            IHmat2 = IHmat;
            IHmat2(4,1) = 0.5*(IHmat(4,2) + IHmat(4,1));
            IH2 = interval(IHmat2(:,1), IHmat2(:,2));
            
            %generate split sets
            R1 = R{iCycle+1} & zonotope(IH1);
            R1 = shrinkIH(R1);
            
            R2 = R{iCycle+1} & zonotope(IH2);
            R2 = shrinkIH(R2);
            
            %choose which split set has to catch up
            if phaseDiff_cycle>0
                RcatchUp = R1;
                Rconst = R2;
            else
                RcatchUp = R2;
                Rconst = R1;
            end
            
            %obtain center difference
            centerDiff = center(R1.Z{2}) - center(R2.Z{2});
            phaseDiff = centerDiff(4);
            
            %number of catchUpCycles
            catchUpCycles = floor(abs(phaseDiff/phaseDiff_cycle));
            
            %loop
            for iCatchUpCycle = 1:catchUpCycles
                %copute reachable set of a single cycle
                RcatchUp = singleCycleReach(A, A_input, U, Uuncertain, t_min, t_max, t_cycle, RcatchUp, options);
            end
            
            %unify
            R{iCycle+1} = enclose(Rconst,RcatchUp);
            %simplify
            R{iCycle+1} = simplifySet(R{iCycle+1}, linMap, A, 1);
            
            %extra cycles
            extraCycles = extraCycles + catchUpCycles;
        end
    end
    
    
%     if mod(iCycle,200)==0
%         dims=[1 2; 2 3; 1 4];
%         for iPlot=1:3
%             figure;
%             plot(R{iCycle+1},dims(iPlot,:));
%         end
%     end
    
    
    if mod(iCycle,100)==0
        disp('next 100 cycles')
        IH
    end

    try
        IH
    catch
    end
    
    iCycle = iCycle+1;
end



%compute phase velocity boundaries
function [vMin,vMax]=phaseVelBound(A,c,vpInt,viInt)

vInt = A(4,:)*[viInt; 0; vpInt; 0] + c(4) + 27;

vMin = inf(vInt);
vMax = sup(vInt);


%compute phase velocity boundaries
function [vMin,vMax]=voltageVelBound(A,vpInt,Ucertain)

%projection matric
P = [1 0 0 0; 0 1 0 0];

absVal = max(abs(inf(vpInt)),abs(sup(vpInt)));

%velocity for v_p1 not considering inputs
vInt_tmp = A(2,:)*[0; -absVal; +absVal; 0];

%velocity for v_i and v_p1 
vIH = interval(P*Ucertain) + interval([0; -vInt_tmp], [0;vInt_tmp]);

vMin = vIH(:,1);
vMax = vIH(:,2);


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


%compute time intervals for charge pump on and off
function [t_on]=timeBound_voltage_conservative(R,vMax,vSat)

%obtain range of v_i and v_p1
v_i_IH = interval([1 0 0 0]*R);
v_p1_IH = interval([0 1 0 0]*R);

v_i_max = v_i_IH(:,2);
v_p1_max = v_p1_IH(:,2);

%t_on
t_i_max = max((0.99*vSat-v_i_max)/vMax(1),0);
t_p1_max = max((vSat-v_p1_max)/vMax(2));
t_on = [t_i_max, t_p1_max];

function factor = inputRange(pllSys, tMax, IH, ind, relInd, options)

%generate corner cases for integral part
V = get(vertices(IH),'V');

for i=1:length(V(1,:))
    %generate initial set
    x0 = zeros(5,1);
    x0(4:5) = [-1 0];
    x0(ind) = V(:,i);
    
    %simulate system
    options.u = options.uLoc{2};
    [obj,t,x] = simulate(pllSys,options,0,tMax,x0,[]);
    
    %find input ranges
    x_min(i) = min(x(:,relInd));
    x_max(i) = max(x(:,relInd));
    
end

xMin = min(x_min);
xMax = max(x_max);

if xMax>0.35
    factor(1) = (1-(xMax-0.35)/(0.51-0.35));
else
    factor(1) = 1;
end

if xMin>0.35
    factor(2) = (1-(xMin-0.35)/(0.51-0.35));
else
    factor(2) = 1;
end




function eAtInt = inputExponential(A,r,options)

%compute Apowers
Apower = powers(A,options);
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


function Apower = powers(A,options)

%initialize 
Apower = cell(1,options.taylorTerms+1);
Apower{1} = A;  
    
%compute powers for each term and sum of these
for i=1:options.taylorTerms
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


function R = simplifySet(R, linMap, A, exactnessFlag)

%do classical reduction
R = reduce(R,'girard',2);

%update options.W
W{1} = linMap;
W{2} = eye(length(A));
% Aenh(2:4,1:4) = A(2:4,1:4);
% Aenh(1,1) = 1;
% Aenh(4,4) = 0.01;
%W{3} = Aenh;
%W{4} = [1 0 0 -1; 0 1 -1 0; 0 1 0 -1; 1 -1 0 0];

%shrink
R = shrink3(R,W,exactnessFlag);

function [E] = expMatRemainder(A, t, options)

%compute auxiliary matrix
tmpMat = 1/factorial(4)*(t^4)*A^4;

%compute E_1 center, generators
E_center = 0.5*tmpMat;
E_gen{1} = 0.5*tmpMat;

for i=5:options.taylorTerms
    tmpMat = 0.5*(t^i)*A^i;
    
    %compute E center, generators
    E_center = E_center + 1/factorial(i)*tmpMat;
    E_gen{end+1} = 1/factorial(i)*tmpMat;
end

%instantiate matrix zonotopes
E_zono = matZonotope(E_center, E_gen);

%compute remainders
Erem = exponentialRemainder(intervalMatrix(0.5*A*t,0.5*A*t),options.taylorTerms);

%final result
E = intervalMatrix(E_zono) + Erem;


function matZinput = inputMatrix(A, A_input, Tmax, tFinal, U, vInt, options)

%initialize
dim = length(A_input);

t_intMat = intervalMatrix(0.5*Tmax,0.5*Tmax); %other time delay: 50e-3

T(1) = infsup(1/2*Tmax,Tmax); %other time delay: 50e-3
T(2) = infsup(1/6*Tmax^2,1/2*Tmax^2); %other time delay: 50e-3
T(3) = infsup(1/8*Tmax^3,1/6*Tmax^3); %other time delay: 50e-3

T_2(1) = infsup(1/2,1); %other time delay: 50e-3
T_2(2) = infsup(1/6,1/2); %other time delay: 50e-3
T_2(3) = infsup(1/8,1/6); %other time delay: 50e-3

%convert to matrix zonotopes
for i=1:length(T)
    gen{1} = rad(T(i));
    t_zonMat{i} = matZonotope(mid(T(i)),gen);
end

%compute matrix zonotopes
E = expMatRemainder(A_input, Tmax, options);

%compute factor
factor = (1/vInt);
matFactor_int = intervalMatrix(mid(factor),rad(factor));
matFactor_zono = matZonotope(matFactor_int);

%zonotope part of input matrix
inputMat_zono = matFactor_zono*(eye(dim) + t_zonMat{1}*A_input...
                + t_zonMat{2}*A_input^2 + t_zonMat{3}*A_input^3);        
inputMat_int = matFactor_int*E*(eye(dim) + 1/factorial(2)*t_intMat*A_input ...
                + 1/factorial(3)*t_intMat^2*A_input^2 + 1/factorial(4)*t_intMat^3*A_input^3 + E);
%----------------------------------------------------------------------

%INPUT-----------------------------------------------------------------


%compute fourth column
U = interval(U);
U_inf = U(:,1);
U_sup = U(:,2);
U_intMat = intervalMatrix(0.5*(U_inf+U_sup), 0.5*(U_sup-U_inf));
fourthCol = (-1)*(intervalMatrix(inputMat_zono) + inputMat_int)*U_intMat;
fourthCol = expm(A*(tFinal-Tmax))*fourthCol; %correction
fourthCol_center = mid(fourthCol.int);
fourthCol_gen = rad(fourthCol.int);

Acenter = zeros(dim);
Adelta = zeros(dim);
Acenter(:,4) = fourthCol_center;
Adelta(:,4) = fourthCol_gen;
matZinput = intervalMatrix(Acenter,Adelta);


function R_new = singleCycleReach(A, A_input, U, Uuncertain, t_min, t_max, t_cycle, Rprev, options)

%compute input solution for t=[0,t_min] 
eAtInt = inputExponential(A_input,t_min,options);
Rinput1 = expm(A*(t_cycle-t_min)) * eAtInt * U;

%compute input solution for t=[t_min,t_max] 
eAtInt = inputExponential(A_input,t_max-t_min,options);
Rinput2 = expm(A*(t_cycle-t_max)) * eAtInt * Uuncertain;

%complete solution
genMat{1} = 0.5*A_input;
B = matZonotope(0.5*A_input,genMat)*(t_max-t_min);
[eBZ,eBI] = expmMixed(B,t_max-t_min,2,6);
eBt = intervalMatrix(eBZ)+eBI;

genMat{1} = 0.5*A;
C = matZonotope(0.5*A,genMat)*(t_max-t_min);
[eCZ,eCI] = expmMixed(C,t_max-t_min,2,6);
eCt = intervalMatrix(eCZ)+eCI;

R_hom = expm(A*(t_cycle-t_max))*eCt*eBt*expm(A_input*t_min)*Rprev;
R_new = R_hom + Rinput1 + Rinput2; 

function D = derivativeSet(A, U)

%define worst case reachable set
Rwc = interval([-10; -10; -10; -pi], [10; 10; 10; pi]);
                
U = interval(U);
                
factor = infsup(0,1);

D = A*interval(Rwc) + factor*interval(U);
D2 = interval(interval(A*zonotope(Rwc))) + factor*interval(U);

%------------- END OF CODE --------------