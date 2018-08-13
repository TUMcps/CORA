function Rguard = partialGuardProjection(obj,hSpace,A,B,loc_old,options)
% partialGuardProjection - computes partial projections of the reachable 
% set onto hyperplanes to separate the part that hits the guard set from
% the part that does not.
%
% Syntax:  
%    Rguard = partialGuardProjection(obj,A,B,t_hit,tmin,tmax,R0,options)
%
% Inputs:
%    obj - location object
%    A - system matrix
%    B - input matrix
%    ...
%    R0 - initial reachable set
%    options - options struct
%
% Outputs:
%    Rguard - intersections with guard set and auxiliary hyperplane
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
% Written:      23-August-2013
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%generate auxiliary hyperplane
c_new = (hSpace.c.'*A).';
d_new = - hSpace.c.'*B*options.u;
h_new = halfspace(c_new, d_new);

%plot hyperplane
plot(h_new);

%continue reachable set computation
%initialize reachable set
contDynamics = get(loc_old,'contDynamics');
[contDynamics,Rnext,options]=initReach(contDynamics,obj,options);

%check if still in front of guard set
intersect = zonoIntersect(h_new, Rnext.ti);

if intersect
    Rlast = R0;
else
    %while reachable set still in front of guard set
    while ~intersect
        %set Rprev
        Rprev = Rnext;

        %compute next reachable set
        [Rnext, options]=reach(contDynamics, Rprev, options);

        %check if still in front of guard set
        intersect = zonoIntersect(h_new, Rnext.ti);
    end
    Rlast = Rprev.tp;
end

%obtain intersection duration
intersectionCounter = 0;
while intersect
    %compute next reachable set
    [Rnext, options]=reach(contDynamics, Rnext, options);

    %check if still in front of guard set
    intersect = zonoIntersect(h_new, Rnext.ti);
    
    %increase counter
    intersectionCounter = intersectionCounter + 1;
end

%obtain t_min, t_max
tmin = 0;
tmax = (intersectionCounter-1)*options.timeStep;

%compute t_hit
%obtain data from location
aux = 1e6*ones(length(c_new),1);
inv = interval(-aux,aux);
tran_old = get(loc_old, 'transition');
reset = get(tran_old{1}, 'reset');
tran{1}=transition(h_new,reset,0,'a','b'); 
%create dummy location
loc = location('dummyLoc',inv,tran,contDynamics);  
eventOptions = odeset('Events',eventFcn(loc));
opt = odeset(eventOptions);
%simulate continuous dynamics
[contDynamics,t] = simulate(contDynamics,options,tmin,tmax,center(Rlast),opt);
t_hit = t(end);


%compute intersection
Rguard = guardProjection(Rlast,h_new,A,B,t_hit,tmin,tmax,contDynamics,options);

%split Rguard with halfspace
Rguard_split = split(Rguard, hSpace);

%region for which guard intersection is ensured for the straight trajectory
%solution
h_region = intersectionRegion(A,B,t_hit,hSpace,options);

%compute backwards solution for Rguard_split{1} and t_max
oldTimeStep = options.timeStep;
options.timeStep = -intersectionCounter*options.timeStep;
options.taylorTerms = 10;
[eAt, eAtInt] = exponentialMatrices(contDynamics,options);
options.timeStep = oldTimeStep;

%R_backwards 
R_backwards = eAt*Rguard_split{1} + eAtInt*options.u;

%hitting time
eventOptions = odeset('Events',eventFcn(loc_old));
opt = odeset(eventOptions);
%simulate continuous dynamics
tmin = 0;
tmax = intersectionCounter*options.timeStep;
[contDynamics,t] = simulate(contDynamics,options,tmin,tmax,center(R_backwards),opt);
t_hit = t(end);

%compute guard intersection
R_pass = guardProjection(R_backwards,hSpace,A,B,t_hit,tmin,tmax,contDynamics,options);


function h_region = intersectionRegion(A,B,t_hit,hSpace,options)

%get data
u = B*options.u;

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

c_region = (hSpace.c.'*(A + Theta)).';
d_region = - hSpace.c.'*Gamma*u;
h_region = halfspace(c_region, d_region);


%------------- END OF CODE --------------