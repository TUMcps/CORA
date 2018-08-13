function R = reachPLL_lock(obj, options)
% reachPLL_lock - computes the reachable set of the PLL when close to
% locking
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
% Written:      03-December-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


%obtain variables
A = options.sys.A;
c = options.sys.c;
dim = length(A);
Ucertain = options.Ucertain;
Ucertain_i = options.Ucertain_i;
Ucertain_p = options.Ucertain_p;
Uuncertain = options.Uuncertain;
R{1} = options.R0; 


%compute phase and voltage velocity boundaries
%set viInt and vpInt
viInt = infsup(-0.5,0.5);
vpInt = infsup(-0.5,0.5);
[vMinPhase,vMaxPhase]=phaseVelBound(A, c, viInt, vpInt);


%compute cycle time
t_cycle = 1/27;

%compute interval hull
IH = interval(R{1});

%compute P_prev
P_prev = polytope(R{1});

iCycle=1;

while ~(IH<=options.Rgoal)

    %PHASE 1---------------------------------------------------------------
    %compute T
    Phi_IH = abs(interval([0 0 0 1]*R{iCycle}))
    Tmax = Phi_IH(:,2)/min(vMinPhase,vMaxPhase);
    
    %compute phase velocity interval
    vInt = infsup(vMinPhase,vMaxPhase);
    
%     %get maximum v_p1 voltage
%     v_IH = interval([0 1 0 0]*R{iCycle});
%     voltagemax = max(abs(v_IH(:,1)), abs(v_IH(:,2)));
%     
%     %compute matZinput
%     if (voltagemax > 0.3) || (iCycle==1)
%         matZinput = inputMatrix(A, Tmax, t_cycle, Ucertain_i, vInt, options);
%     else
        matZinput = inputMatrix(A, Tmax, t_cycle, Ucertain, vInt, options);
%    end

    
    %complete solution
    Rinit_new = (expm(A*t_cycle)+matZinput)*R{iCycle};
    
    %intersect with invariant
    R{iCycle+1} = Rinit_new & zonotope(options.Rinv);
    R{iCycle+1} = shrinkIH(R{iCycle+1});
    R{iCycle+1} = reduce(R{iCycle+1},'girard',60);
    
    %if mod(iCycle,200)==0
    %if mod(iCycle,100)==0
    if mod(iCycle,2000)==0 
        %obtain linear map
        t = iCycle*t_cycle;
        linMap = expm(options.sys.A*t) + 0.01*eye(dim);
        
        %simplify
        R{iCycle+1} = simplifySet(R{iCycle+1}, linMap, A, 1);
    end
    
    
%     if mod(iCycle,200)==0
%         dims=[1 2; 2 3; 1 4];
%         for iPlot=1:3
%             figure;
%             plot(R{iCycle+1},dims(iPlot,:));
%         end
%     end
    
    IH=interval(R{iCycle+1});
    
%     if mod(iCycle,50)==0
%         save lockModeAll
%     end

    iCycle = iCycle+1;
end
R


%compute phase velocity boundaries
function [vMin,vMax]=phaseVelBound(A,c,vpInt,viInt)

vInt = A(4,:)*[viInt; 0; vpInt; 0] + c(4);

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
vIH = interval(P*Ucertain) + interval([0; -vInt_tmp],[0; vInt_tmp]);

vMin = vIH(:,1);
vMax = vIH(:,2);


%compute time intervals for charge pump on and off
function [t_on,t_total]=timeBound_phase(R,vMin,vMax,phaseMargin)

%obtain range of Phi_v
Phi_IH = interval([0 0 0 1]*R);
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


function R = simplifySet(R, linMap, A, exactnessFlag)

%do classical reduction
R = reduce(R,'girard',2);

%update options.W
W{1} = eye(length(A));
W{2} = linMap;
% Aenh(2:4,1:4) = A(2:4,1:4);
% Aenh(1,1) = 1;
% Aenh(4,4) = 0.01;
%W{3} = Aenh;
%W{4} = [1 0 0 -1; 0 1 -1 0; 0 1 0 -1; 1 -1 0 0];

%shrink
R = shrink3(R,W,exactnessFlag);

function [E] = expMatRemainder(A, t, options)

%compute auxiliary matrix
tmpMat = 0.5*(t^2)*A^2;

%compute E_1 center, generators
E_center = 1/factorial(2)*tmpMat;
E_gen{1} = 1/factorial(2)*tmpMat;

for i=3:options.taylorTerms
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


function matZinput = inputMatrix(A, Tmax, tFinal, U, vInt, options)

%initialize
dim = length(A);

T_1 = infsup(0,Tmax); %other time delay: 50e-3
T_2 = infsup(0.5*Tmax,Tmax); %other time delay: 50e-3

%get center and delta of [0,T]
t_center_1 = mid(T_1);
t_delta_1 = rad(T_1);
t_intMat_1 = intervalMatrix(t_center_1,t_delta_1);

%get center and delta of [0.5T, T]
t_center_2 = mid(T_2);
t_delta_2{1} = rad(T_2);
t_zonMat_2 = matZonotope(t_center_2,t_delta_2);

%compute matrix zonotopes
E = expMatRemainder(A, Tmax, options);

%compute factor
factor = (1/vInt);
matFactor_int = intervalMatrix(mid(factor),rad(factor));
matFactor_zono = matZonotope(matFactor_int);

%zonotope part of input matrix
inputMat_zono = matFactor_zono*(eye(dim) + t_zonMat_2*A);
inputMat_int = matFactor_int*(E + E*0.5*A*t_intMat_1 + E*E);
%----------------------------------------------------------------------

%INPUT-----------------------------------------------------------------


u_IH = get(interval(U) + (-options.sys.c),'intervals');
u_IH = -u_IH; %negative feedback
u_center = 0.5*(u_IH(:,1)+u_IH(:,2));
u_delta{1} =  0.5*(u_IH(:,2)-u_IH(:,1));
u_int = intervalMatrix(u_center, u_delta{1});
u_zono = matZonotope(u_center, u_delta);

%compute fourth column
fourthCol = intervalMatrix(inputMat_zono*u_zono) + inputMat_int*u_int;
fourthCol = expm(A*(tFinal-Tmax))*fourthCol; %correction
fourthCol_center = mid(fourthCol.int);
fourthCol_gen = rad(fourthCol.int);

Acenter = zeros(dim);
Adelta = zeros(dim);
Acenter(:,4) = fourthCol_center;
Adelta(:,4) = fourthCol_gen;
matZinput = intervalMatrix(Acenter,Adelta);


%------------- END OF CODE --------------