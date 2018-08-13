function Rguard = guardProjection_constFlow(obj,halfspace,A,B,t_hit,tmin,tmax,options)
% hyperplaneMap - computes the reachable set of the system within a 
% location, detects the guard set that is hit and computes the new set on
% the hyperplane and the subsequent mapping.
%
% Syntax:  
%    [TP,R,activeGuards,Rjump,Rcont] =
%    hyperplaneMap(obj,tStart,R0,options)
%
% Inputs:
%    obj - location object
%    tStart - start time
%    R0 - initial reachable set
%    options - options struct
%
% Outputs:
%    TP - time point struct; e.g. contains time vector of minimum times for reaching guard sets
%    R - cell array of reachable sets
%    activeGuards - active guards
%    Rjump - reachable set after jump according to the reset map
%    Rcont - reachable set due to continuous evolution without guard or
%    invariant inmtersection
%
% Example: 
%
% Other m-files required: initReach, reach, potOut, potInt, guardIntersect,
% reset
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      29-August-2013
% Last update:  25-July-2016 (intervalhull replaced by interval)---
% Last revision:---

%------------- BEGIN CODE --------------

%get data
x0 = center(obj);
u = B*options.u;
n = halfspace.c;
d = halfspace.d;

%compute f
dim = length(A);
I = eye(dim);
%Theta
Theta = A^2*t_hit/2;
for i=3:options.taylorTerms
    Theta = Theta + A^i*t_hit^(i-1)/factorial(i);
end
%Gamma
Gamma = I;
for i=1:(options.taylorTerms-1)
    Gamma = Gamma + A^i*t_hit^i/factorial(i+1);
end
f = (A+Theta)*x0 + Gamma*u;

%use reduced zonotope 
Rred = reduce(obj,'girard',options.errorOrder);

%projection on halfspace
Rhit_const = (I - f*n.'/(n.'*f))*Rred + f*d/(n.'*f);

%for test
Zk = Rhit_const + zonotope([zeros(dim,1), 0.01*eye(dim)]);

%compute mapping error
delta_Rhit = mappingError(A,f,u,obj,tmin,tmax,t_hit,options.taylorTerms,halfspace.c);
%delta_Rhit = mappingError_split(A,f,u,obj,tmin,tmax,t_hit,options.taylorTerms,Lambda_int,halfspace.c,options.maxTimeStep);

%compute set on hyperplane
Rguard = Rhit_const + delta_Rhit;






%NOT YET UPDATED!
function delta_Rhit = mappingError_split(A,f,u,R0,tmin,tmax,tc,order,Lambda_int,n,tmaxStep)

%obtain center and delta set
x0 = center(R0);

%compute time span
timeSpan = tmax - tmin;

%obtain necessary subintervals
subIntervals = ceil(timeSpan/tmaxStep);

%obtain time intervals
delta_t = timeSpan/subIntervals;

for iTimeInterval = 1:subIntervals
    %obtain new lower and upper bounds
    tmin_new = tmin + (iTimeInterval-1)*delta_t;
    tmax_new = tmin + iTimeInterval*delta_t;
    
    %compute powers of tmin and tmax
    tZono = timePowers(tmin_new,tmax_new,order);
 
    %auxiliary set
    Raux_1 = zonotope(0*x0);
    Raux_2 = zonotope(0*x0);
    for i=2:order
        Raux_1 = Raux_1 + A^i/factorial(i)*(tZono{i-1}*R0 + (-tc^(i-1)*x0));
    end
    for i=1:order
        Raux_2 = Raux_2 + A^i/factorial(i+1)*(tZono{i}+(-tc^i))*u;
    end
    
    if iTimeInterval==1
        R_union = interval(Raux_1 + Raux_2);
    else
        R_union = R_union | interval(Raux_1 + Raux_2);
    end
end
% convert to zonotope
R_union = zonotope(R_union);

%compute time powers of full time interval
tZono = timePowers(tmin,tmax,order);

%compute remainder matrix
A_abs = abs(A);
Apower_abs = eye(length(A_abs));
M = eye(length(A_abs));

for i=1:order
    Apower_abs = Apower_abs*A_abs;
    M = M + Apower_abs*tmax^i/factorial(i);
end 
  
%determine error due to finite Taylor series
W = expm(A_abs*tmax) - M;
%compute absolute value of W for numerical stability
W = abs(W);
E = interval(-W,W);

%final additional set
G = tZono{1}*(R_union) + E*R0 + E*tmax*u;

%error due to uncertain hitting time
delta_t = interval((-n')*G)/Lambda_int;
delta_t_mid = mid(delta_t);
delta_t_rad{1} = rad(delta_t);
delta_t_zono = matZonotope(delta_t_mid,delta_t_rad);

%final uncertainty
delta_Rhit = G + delta_t_zono*(A*R0 + f);



function delta_Rhit = mappingError(A,f,u,R0,tmin,tmax,tc,order,n)

%compute powers of tmin and tmax
tZono = timePowers(tmin,tmax,order);

%obtain center and delta set
x0 = center(R0);
%dR = R0 + (-x0);

%auxiliary set: new
Raux_1 = A*(R0 + (-x0)); %<-- changed for constant flow!
Raux_2 = zonotope(0*x0);


for i=2:order 
    Raux_1 = Raux_1 + A^i/factorial(i)*(tZono{i-1}*R0 + (-tc^(i-1)*x0));
end
for i=1:order
    Raux_2 = Raux_2 + A^i/factorial(i+1)*(tZono{i}+(-tc^i))*u;
end

%compute remainder matrix
A_abs = abs(A);
Apower_abs = eye(length(A_abs));
M = eye(length(A_abs));

for i=1:order
    Apower_abs = Apower_abs*A_abs;
    M = M + Apower_abs*tmax^i/factorial(i);
end 
  
%determine error due to finite Taylor series
W = expm(A_abs*tmax) - M;
%compute absolute value of W for numerical stability
W = abs(W);
E = interval(-W,W);

%final additional set
G = tZono{1}*(Raux_1 + Raux_2) + E*R0 + E*tmax*u;

%final uncertainty
delta_Rhit = G + (-f*n.'/(n.'*f)*G);


function tZono = timePowers(tmin,tmax,order)


%first order
tminPow(1) = tmin;
tmaxPow(1) = tmax;

tradPow(1) = 0.5*(tmaxPow(1) - tminPow(1));
tmidPow(1) = 0.5*(tmaxPow(1) + tminPow(1));

tRad{1} = tradPow(1);
tZono{1} = matZonotope(tmidPow(1),tRad);

%higher orders
for i=1:order
    tminPow(i+1) = tminPow(i)*tmin;
    tmaxPow(i+1) = tmaxPow(i)*tmax;
    
    tradPow(i+1) = 0.5*(tmaxPow(i+1) - tminPow(i+1));
    tmidPow(i+1) = 0.5*(tmaxPow(i+1) + tminPow(i+1));
    
    tRad{1} = tradPow(i+1);
    tZono{i+1} = matZonotope(tmidPow(i+1),tRad);
end




%------------- END OF CODE --------------