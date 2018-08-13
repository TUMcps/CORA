function [Rfinal] = reachPLL_cycles(Rinit, cycles, options)
% reachPLL_cycles - computes the reachable set of the PLL for a fixed
% number of cycles
%
% Syntax:  
%    [R] = reachPLL_cycles(options)
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
% Written:      02-December-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


%obtain variables
A = options.sys.A;
c= options.sys.c;
dim = length(A);
Ucertain = options.Ucertain;
Ucertain_i = options.Ucertain_i;
Ucertain_p = options.Ucertain_p;
Uuncertain = options.Uuncertain;
guard = options.guard;
R = cell(1,cycles+1);
R{1} = Rinit;

%extract minimum phase margin
IHmargin = interval([0 0 0 1 0]*options.Rgoal{1});
phaseMargin = IHmargin(:,2);

%instantiate linear system
linSys = linearSys('linearDynamics',A,eye(dim)); 

%compute phase and voltage velocity boundaries
%set viInt and vpInt
viInt = infsup(-1,1);
vpInt = infsup(-1,1);
vSat = 0.7; %<-- change vSat here!
[vMinPhase,vMaxPhase]=phaseVelBound(A, c, viInt, vpInt);
[vMinVoltage,vMaxVoltage] = voltageVelBound(A,vpInt,Ucertain);

%compute cycle time
t_cycle = 1/options.sys.c(5);

%compute interval hull
IH = interval(R{1});

for iCycle=1:cycles

    %compute minimum and maximum times for locations
    [t_on_phase,t_total]=timeBound_phase(R{iCycle},vMinPhase,vMaxPhase,phaseMargin);
    try
    [t_on_voltage]=timeBound_voltage(R{iCycle},vMaxVoltage,vSat);
    catch
        disp('error');
    end
    %choose minimum time for t_on
    t_on_min_vi = min(inf(t_on_phase), t_on_voltage(1));
    t_on_min_vp1 = min(inf(t_on_phase), t_on_voltage(2));
    t_on_max = sup(t_on_phase);

    %compute time increments without considering saturation
    t1 = min(t_on_min_vi, t_on_min_vp1);
    t2 = max(t_on_min_vi, t_on_min_vp1);
    t3 = t_on_max;
    t4 = t_cycle;
    deltaT(1) = t1;
    deltaT(2) = t2-t1;
    deltaT(3) = t3-t2;
    deltaT(4) = t4-t3;
    
    %check correctness
    if (deltaT(1) < 0) || (deltaT(1)+deltaT(2)+deltaT(3)+deltaT(4)~=t_cycle)
        disp('time error!!');
    end

    %compute mapping matrix when in location 2
    M = cell(1,4);
    eAtInt = cell(1,4);
    for i=1:4
        M{i} = expm(options.sys.A*deltaT(i));
        eAtInt{i} = inputExponential(A,deltaT(i),options);
    end
    
    %determine uncertain inputs
    U = cell(1,4);
    U{1} = Ucertain;
    if t_on_voltage(1) > t_on_voltage(2)
        U{2} = Ucertain_i;
    else
        U{2} = Ucertain_p;
    end
    U{3} = Uuncertain;
    U{4} = zonotope(c);
    

    %check if reachable set is still in non-saturation
    intersectIH = IH & options.RnonSat;
    interval = get(intersectIH,'intervals');
    if isempty(interval)
        disp('saturation assumption violated!!!------------------')
    end

    %compute reachable set from t=0 to t=t_cycle when using location 2

    % t = [0, deltaT_1]
    Rinput = cell(1,4);
    Rtmp = cell(1,5);
    Rinput{1} = eAtInt{1}*U{1};
    Rtmp{1} = M{1}*R{iCycle} + Rinput{1};
    
    for k=2:4
        Rinput{k} = eAtInt{k}*U{k};
        Rtmp{k} = M{k}*Rtmp{k-1} + Rinput{k};
    end
    
    %reset
    Rtmp{5} = Rtmp{4} + [0; 0; 0; -1; -1];

    %reduce zonotope
    R{iCycle+1}=reduce(Rtmp{end},'girard',options.zonotopeOrder);
    
    %compute interval hull
    IH = interval(R{iCycle+1});
    
    %check enclosure in invariant
    if ~(IH <= options.Rinv)
        %intersect
        R{iCycle+1} = R{iCycle+1} & zonotope(options.Rinv);
        R{iCycle+1} = shrinkIH(R{iCycle+1});
    end
    
    
    %compute reachable set intersecting guard set 2
    if inf(t_total)<t_cycle && options.overflow==1
        %compute reachable set at time t_total_min
        deltaT_other = inf(t_total) - t3;
        M_other = expm(options.sys.A*deltaT_other);
        eAtInt_other = inputExponential(A,deltaT_other,options);
        
        %compute R_init
        Rinput_other = eAtInt_other*U{4};
        Rinit = M_other*Rtmp{3} + Rinput_other;
        
        %set options
        options.timeStep = sup(t_total) - inf(t_total);
        options.uTrans = center(U{4});
        options.U = U{4} + (-options.uTrans);
        
        %compute initial set
        [obj,Rfirst,options] = initReach(linSys,Rinit,options);
        
        %plot results
%         plot(Rinit, [4 5]);
%         plot(Rfirst.ti, [4 5]);
%         plot(guard, [4 5]);
    end
end
%store final set
Rfinal = R{iCycle+1};



%compute phase velocity boundaries
function [vMin,vMax]=phaseVelBound(A,c,vpInt,viInt)

vInt = A(4,:)*[viInt; 0; vpInt; 0; 0] + c(4);

vMin = inf(vInt);
vMax = sup(vInt);


%compute phase velocity boundaries
function [vMin,vMax]=voltageVelBound(A,vpInt,Ucertain)

%projection matric
P = [1 0 0 0 0; 0 1 0 0 0];

absVal = max(abs(inf(vpInt)),abs(sup(vpInt)));

%velocity for v_p1 not considering inputs
vInt_tmp = A(2,:)*[0; -absVal; +absVal; 0; 0];

%velocity for v_i and v_p1 
vIH = interval(P*Ucertain) + interval([0;-vInt_tmp], [0; vInt_tmp]);

vMin = vIH(:,1);
vMax = vIH(:,2);


%compute time intervals for charge pump on and off
function [t_on,t_total]=timeBound_phase(R,vMin,vMax,phaseMargin)

%obtain range of Phi_v
Phi_IH = interval([0 0 0 1 0]*R);
PhiMin = min(Phi_IH(:,1),phaseMargin);
PhiMax = min(Phi_IH(:,2),phaseMargin);

%t_on
t_on_min = -PhiMax/vMax;
t_on_max = -PhiMin/vMin;
t_on = infsup(t_on_min,t_on_max);


%t_total
t_total_min = (1-PhiMax)/vMax;
t_total_max = (1-PhiMin)/vMin;
t_total = infsup(t_total_min,t_total_max);


%compute time intervals for charge pump on and off
function [t_on]=timeBound_voltage(R,vMax,vSat)

%obtain range of v_i and v_p1
v_i_IH = interval([1 0 0 0 0]*R);
v_p1_IH = interval([0 1 0 0 0]*R);

v_i_min = v_i_IH(:,1);
v_p1_min = v_p1_IH(:,1);

%t_on
t_i_max = (vSat-v_i_min)/vMax(1);
t_p1_max = (vSat-v_p1_min)/vMax(2);
t_on = [t_i_max, t_p1_max];



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


%------------- END OF CODE --------------