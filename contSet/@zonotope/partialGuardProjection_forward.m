function Rguard = partialGuardProjection_forward(obj,hSpace,A,B,loc_old,options)
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
% Last update:  27-August-2013
% Last revision:---

%------------- BEGIN CODE --------------

%generate auxiliary hyperplane
c_new = (hSpace.c.'*A).';
d_new = - hSpace.c.'*B*options.u;
h_new = halfspace(c_new, d_new);

%plot hyperplane
plot(h_new);

%determine switching time interval
[tmin,tmax,t_hit,t_total,Rlast] = switchingTimeInterval(obj,loc_old,h_new,options);

%compute intersection
contDynamics = get(loc_old,'contDynamics');
Rguard = guardProjection(Rlast,h_new,A,B,t_hit,tmin,tmax,loc_old,0,options);

%compute rotated halfspace for better accuracy 
%find rotation Point
rotPoint = commonPoint(hSpace, h_new);
%obtain new direction
newDir = hSpace.c - hSpace.c.'*h_new.c/norm(h_new.c)^2*h_new.c;
%perform rotation
hRot = rotate(hSpace, newDir, rotPoint);

%split Rguard with halfspace
Rguard_split_old = split(Rguard, hSpace);
Rguard_split = split(Rguard, hRot);

%region for which guard intersection is ensured for the straight trajectory
%solution
h_region = intersectionRegion(A,B,t_hit,hSpace,options);

%compute backwards solution for Rguard_split{1} and t_max
oldTimeStep = options.timeStep;
options.timeStep = -t_total;
options.taylorTerms = 10;
[eAt, eAtInt] = exponentialMatrices(contDynamics,options);
options.timeStep = oldTimeStep;

%R_backwards 
R_backwards = eAt*Rguard_split{1} + eAtInt*options.u;

%hitting time
[tmin,tmax,t_hit,t_total,Rlast,x_hit] = switchingTimeInterval(R_backwards,loc_old,hSpace,options);

%obtain intermediate switching surfaces
%obtain necessary sub-halfspaces
%nrOfHalfspaces = ceil(tmax/options.maxTimeStep);
nrOfHalfspaces = 2;
x_first = center(R_backwards);
x_last = x_hit;
dir_first = x_first - x_last;
dir_first = dir_first/norm(dir_first);
dir_last = hSpace.c;
dir_last = dir_last/norm(dir_last);
hIntermediate = intermediateHalfspaces(x_first, x_last, dir_first, dir_last, hSpace, nrOfHalfspaces);

% %move reachable set further back
% oldTimeStep = options.timeStep;
% options.timeStep = -0.1*t_total;
% options.taylorTerms = 10;
% [eAt, eAtInt] = exponentialMatrices(contDynamics,options);
% options.timeStep = oldTimeStep;
% 
% %R_backwards 
% R_backwards = eAt*R_backwards + eAtInt*options.u;


%compute guard intersections
R_pass = R_backwards;
for i = 1:length(hIntermediate)
    
    [tmin,tmax,t_hit,t_total,R_pass] = switchingTimeIntervalAdvanced(R_pass,loc_old,hIntermediate{i},h_new,A,B,options);

    %compute guard projection
    R_pass = guardProjection(R_pass,hIntermediate{i},A,B,t_hit,tmin,tmax,loc_old,0,options);

    %project result
    R_pass = project(hIntermediate{i}, R_pass);

    %reduce order
    R_pass = reduce(R_pass,'girard',options.zonotopeOrder);
    
    %plot
    plot(hIntermediate{i});
    Rred = reduce(R_pass, 'girard', 10);
    Rred = Rred + zonotope([zeros(2,1), 1e-4*eye(2)]);
    plot(Rred, [1 2], 'r');
end



function [tmin,tmax,t_hit,t_total,Rlast] = switchingTimeIntervalAdvanced(R0,loc,h,h_sep,A,B,options)

%compute flow in the center of R0
f = A*center(R0) + B*options.u;

%compute movement after a step size of options.maxTimeStep
factor = 0.1;
delta_x = factor*f*options.maxTimeStep;

%init 1-dimensional partitioning
intersect = 1;
%h_sep_trans = halfspace(-h_sep.c,h_sep.d); %reverse direction of normal vector
h_sep_trans = h_sep; 
R_rem = R0;
iStep = 1;

%as long as the halfspace has not completely moved over R0
while intersect
    %translate halfspace
    h_sep_trans = h_sep_trans + (-delta_x);
    
    %check if halfspace still intersects R0
    intersect = zonoPartiallyIn(h_sep_trans, R0);
    
    %compute rotated halfspace for better accuracy 
    %find rotation Point
    rotPoint = commonPoint(h_sep_trans, h);
    %obtain new direction
    newDir = h_sep_trans.c - h_sep_trans.c.'*h.c/norm(h.c)^2*h.c;
    %perform rotation
    hRot = rotate(h_sep_trans, newDir, rotPoint);
    
    %compute split
    %R_split = split(R_rem, h_sep_trans);
    R_split = split(R_rem, hRot);
    
    if iscell(R_split)
        R_rem = R_split{1};
        R_partial{iStep} = R_split{2};

        %update iStep
        iStep = iStep + 1;
    end
end
%last remaining set becomes a partial set
R_partial{iStep} = R_rem;

%compute t_min and t_max for each partial reachable set
tmin = inf;
tmax = -inf;
for iSet = 1:length(R_partial)
    [tmin_1,tmax_1,t_total,Rlast] = intersectionTime(R_partial{iSet},loc,h,options);
    [tmin_2,t_dummy,tmax_2] = intersectionTime(Rlast,loc,h_sep,options);
    
    tmin_iSet = min(tmin_1, tmax_2);
    tmax_iSet = min(tmax_1, tmax_2);
    
    tmin = min(tmin, tmin_iSet);
    tmax = max(tmax, tmax_iSet);
end

%compute Rlast and t_total
[t_dummy,t_dummy,t_total,Rlast] = intersectionTime(R0,loc,h,options);


%get contDynamics
contDynamics = get(loc,'contDynamics');

%compute t_hit
%obtain data from location
aux = 1e6*ones(length(options.x0),1);
inv = interval(-aux, aux);
tran_old = get(loc, 'transition');
reset = get(tran_old{1}, 'reset');
tran{1}=transition(h,reset,0,'a','b'); 
%create dummy location
loc = location('dummyLoc',inv,tran,contDynamics);  
eventOptions = odeset('Events',eventFcn(loc));
opt = odeset(eventOptions);
%simulate continuous dynamics
[contDynamics,t,x] = simulate(contDynamics,options,tmin,tmax,center(R0),opt);
t_hit = t(end) - tmin;

%reverse dynamics
if t(end) == tmax
    A = get(contDynamics,'A'); %get system matrix
    B = get(contDynamics,'B'); %get input matrix
    contDynamics_rev = linearSys('reverse_Sys',-A,-B); %initialize linear interval system;
    %simulate continuous dynamics
    [contDynamics,t,x] = simulate(contDynamics_rev,options,tmin,tmax,center(R0),opt);
    t_hit = -(t(end)-tmin);
    
    % initial state exactly on switching surface
    if t(end) == tmax
        t_hit = 0;
    end
end





function [tmin,tmax,t_hit,t_total,Rlast,x_hit] = switchingTimeInterval(obj,loc_old,h,options)

%get intersection times
[tmin,tmax,t_total,Rlast] = intersectionTime(obj,loc_old,h,options);

%get contDynamics
contDynamics = get(loc_old,'contDynamics');

%compute t_hit
%obtain data from location
aux = 1e6*ones(length(options.x0),1);
inv = interval(-aux,aux);
tran_old = get(loc_old, 'transition');
reset = get(tran_old{1}, 'reset');
tran{1}=transition(h,reset,0,'a','b'); 
%create dummy location
loc = location('dummyLoc',inv,tran,contDynamics);  
eventOptions = odeset('Events',eventFcn(loc));
opt = odeset(eventOptions);
%simulate continuous dynamics
[contDynamics,t,x] = simulate(contDynamics,options,tmin,tmax,center(Rlast),opt);
t_hit = t(end) - tmin;
x_hit = x(end,:).';



function [tmin,tmax,t_total,Rlast] = intersectionTime(obj,loc_old,h,options)

%initialize reachable set
contDynamics = get(loc_old,'contDynamics');

%positive direction
[tmin_pos,tmax_pos,t_total_pos,Rlast_pos] = intersectionTime_oneDirection(obj,contDynamics,h,options);

%reverse dynamics
A = get(contDynamics,'A'); %get system matrix
B = get(contDynamics,'B'); %get input matrix
contDynamics_rev = linearSys('reverse_Sys',-A,-B); %initialize linear interval system;

%positive direction
[tmin_neg,tmax_neg,t_total_neg,Rlast_neg] = intersectionTime_oneDirection(obj,contDynamics_rev,h,options);

%if negative direction is not empty 
if ~isempty(tmin_neg)
    if ~isempty(tmin_pos)
        tmin = -max(tmax_neg);
        tmax = max(tmax_pos);
        t_total = [];
        Rlast = Rlast_pos;
    else
        tmin = -tmax_neg;
        tmax = -tmin_neg;
        t_total = [];
        Rlast = Rlast_neg;
    end
else
    tmin = tmin_pos;
    tmax = tmax_pos;
    t_total = t_total_pos;
    Rlast = Rlast_pos;
end




function [tmin,tmax,t_total,Rlast] = intersectionTime_oneDirection(obj,contDynamics,h,options)

%refine
options.timeStep = 0.1*options.timeStep;

%initialize reachable set
[contDynamics,Rnext,options]=initReach(contDynamics,obj,options);

%is distance incearsing?
delta_d = abs(h.c.'*center(Rnext.tp)-h.d) - abs(h.c.'*center(obj)-h.d);
if delta_d > 0
    dist_incr = 1; 
else
    dist_incr = 0; 
end

%check if still in front of guard set
intersect = zonoIntersect(h, Rnext.ti);

if ~dist_incr || intersect

    noIntersectionCounter = 0;
    if intersect
        Rlast = obj;
    else
        %while reachable set still in front of guard set
        while ~intersect
            %set Rprev
            Rprev = Rnext;

            %compute next reachable set
            [Rnext, options]=reach(contDynamics, Rprev, options);

            %check if still in front of guard set
            intersect = zonoIntersect(h, Rnext.ti);

            %increase counter
            noIntersectionCounter = noIntersectionCounter + 1;
        end
        Rlast = Rprev.tp;
    end

    %obtain intersection duration
    intersectionCounter = 0;
    while intersect
        %compute next reachable set
        [Rnext, options]=reach(contDynamics, Rnext, options);

        %check if still in front of guard set
        intersect = zonoIntersect(h, Rnext.ti);

        %increase counter
        intersectionCounter = intersectionCounter + 1;
    end

    %obtain t_min, t_max
    tmin = 0;
    tmax = intersectionCounter*options.timeStep;
    t_total = (intersectionCounter + noIntersectionCounter - 1)*options.timeStep;
    
else
    tmin = [];
    tmax = [];
    t_total = [];
    Rlast = [];
end



function hIntermediate = intermediateHalfspaces(x_first, x_last, dir_first, dir_last, h, nrOfHalfspaces)

%delta position and orientation
delta_pos = (x_last - x_first)/2;
delta_dir = (dir_last - dir_first)/2;

%current halfspace
curr_h = h;

curr_pos = x_first;
curr_dir = dir_first;

%shift halfspace
curr_h = h + (-2*delta_pos);
%rotate halfspace
curr_h = rotate(curr_h, curr_dir, curr_pos);

for i = 1:nrOfHalfspaces
    
    %update position and direction
    curr_pos = curr_pos + delta_pos;
    curr_dir = curr_dir + delta_dir;
    
    %shift halfspace
    curr_h = curr_h + delta_pos;
    
    %rotate halfspace
    curr_h = rotate(curr_h, curr_dir, curr_pos);
    
    %save halfspace
    hIntermediate{i} = curr_h;
    
    %update delta_pos and delta_dir
    delta_pos = delta_pos/2;
    delta_dir = delta_dir/2;
end
%save actual halfspace
hIntermediate{nrOfHalfspaces + 1} = h;



function hIntermediate = intermediateHalfspaces_old(x_first, x_last, dir_first, dir_last, h, nrOfHalfspaces)

%delta position and orientation
delta_pos = (x_first - x_last)/nrOfHalfspaces;
delta_dir = (dir_first - dir_last)/nrOfHalfspaces;

%current halfspace
curr_h = h;
curr_pos = x_last;
curr_dir = dir_last;
for i = 1:nrOfHalfspaces
    
    %update position and direction
    curr_pos = curr_pos + delta_pos;
    curr_dir = curr_dir + delta_dir;
    
    %shift halfspace
    curr_h = curr_h + delta_pos;
    
    %rotate halfspace
    curr_h = rotate(curr_h, curr_dir, curr_pos);
    
    %save halfspace
    hIntermediate{nrOfHalfspaces - i + 1} = curr_h;
end
%save actual halfspace
hIntermediate{nrOfHalfspaces + 1} = h;




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