function [R_guard_noInt, R_new] = partialGuardProjection(obj,hSpace,A,B,loc_old,options)
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
c_new = - (hSpace.c.'*A).';
d_new = hSpace.c.'*B*options.u;
h_new = halfspace(c_new, d_new);

if isfield(options,'debug') && options.debug
    %plot hyperplane
    plot(h_new);
end

%compute intersection
contDynamics = get(loc_old,'contDynamics');
Rguard = guardProjection(obj,h_new,A,B,loc_old,0,options);

%split Rguard with halfspace
Rguard_split = split(Rguard{1}, hSpace);

% %region for which guard intersection is ensured for the straight trajectory
% %solution
% h_region = intersectionRegion(A,B,t_hit,hSpace,options);

%R_backwards 
R_backwards = Rguard_split{1};
R_guard_noInt = Rguard_split{2};

if all(options.projectedDimensions == [1 2])
    %plot
    if options.debug
        dim = length(center(R_backwards));
        Rred = reduce(R_backwards, 'girard', 10);
        Rred = Rred + zonotope([zeros(dim,1), 1e-4*eye(dim)]);
        plot(Rred, [1 2], 'r');
    end
end

%reverse dynamics
A = get(contDynamics,'A'); %get system matrix
B = get(contDynamics,'B'); %get input matrix
A_rev = -A;
B_rev = -B;
contDynamics_rev = linearSys('reverse_Sys',-A,-B); %initialize linear interval system;

%create loc_rev with reverse dynamics
%obtain data from location
aux = 1e6*ones(length(options.x0),1);
inv = interval(-aux, aux);
tran_old = get(loc_old, 'transition');
reset = get(tran_old{1}, 'reset');
tran{1}=transition(hSpace,reset,0,'a','b'); 
%create dummy location
loc_rev = location('reverseLoc',0,inv,tran,contDynamics_rev);  

%obtain intermediate switching surfaces
%obtain necessary sub-halfspaces
%nrOfHalfspaces = ceil(tmax/options.maxTimeStep);
nrOfHalfspaces = 1;
hIntermediate{1} = hSpace;
%hIntermediate = intermediateHalfspaces_linear(h_new, hSpace, rotPoint, nrOfHalfspaces);
%hIntermediate = intermediateHalfspaces(h_new, hSpace, rotPoint, nrOfHalfspaces);


%compute guard intersections
R_new{1} = R_backwards;
for i = 1:length(hIntermediate)
    
    %set counter
    counter = 1;
    while counter <= length(R_new)
    
        %if split is not required
        splitFlag = 1;
        while splitFlag

            %compute guard projection
            [R_guard, dummy, splitFlag] = guardProjection(R_new{counter},hIntermediate{i},A_rev,B_rev,loc_rev,0,options);
            
            if ~isempty(R_guard{1}) && zonoPartiallyIn(hSpace, R_new{counter})

                if splitFlag
                    if i==1
                        R_aux = split(R_new{counter},h_new.c,hIntermediate{1}.c);
                        %project
                        R_aux{1} = project(h_new,R_aux{1});
                        R_aux{2} = project(h_new,R_aux{2});
                    else
                        R_aux = split(R_new{counter},hIntermediate{i-1}.c,h_new.c);
                        %project
                        R_aux{1} = project(hIntermediate{i-1},R_aux{1});
                        R_aux{2} = project(hIntermediate{i-1},R_aux{2});
                    end

                    %delete current R_new
                    R_new(counter) = [];
                    %determine number of sets
                    nrOfSets = length(R_new);
                    %insert split sets
                    R_new{nrOfSets + 1} = R_aux{1};
                    R_new{nrOfSets + 2} = R_aux{2};
                else
                    R_new{counter} = R_guard{1};
                    counter = counter + 1;

                    if options.debug
                        %plot
                        Rred = project(R_new{counter - 1}, [1 2]);
                        Rred = reduce(Rred, 'girard', 10);
                        Rred = Rred + zonotope([zeros(2,1), 1e-4*eye(2)]);
                        plot(Rred, [1 2], 'g');
                    end
                end
                
            else
                %delete current R_new
                R_new(counter) = [];
            end

        end
    end
    


    
%     %plot
%     plot(hIntermediate{i});
%     for i = 1:length(R_new)
%         Rred = reduce(R_new{i}, 'girard', 10);
%         Rred = Rred + zonotope([zeros(2,1), 1e-4*eye(2)]);
%         if mod(i,2) == 0
%             plot(Rred, [1 2], 'r');
%         else
%             plot(Rred, [1 2], 'g');
%         end
%     end
end

%TODO: chop off "left side"



function [tmin,tmax,t_hit,t_total,RbeforeInt] = switchingTimeIntervalAdvanced(R0,loc,h,h_sep,A,B,options)

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

%get contDynamics
contDynamics = get(loc,'contDynamics');

%compute t_min and t_max for each partial reachable set
tmin = inf;
tmax = -inf;
for iSet = 1:length(R_partial)
    %[tmin_1,tmax_1,t_total,RbeforeInt] = intersectionTime(R_partial{iSet},loc,h,options);
    [tmin_1,tmax_1,t_total,RbeforeInt] = intersectionTime_oneDirection(h,R_partial{iSet},contDynamics,options);
    
    tmin = min(tmin, tmin_1);
    tmax = max(tmax, tmax_1);
end

%compute RbeforeInt and t_total
[t_dummy,t_dummy,t_total,RbeforeInt] = intersectionTime(R0,loc,h,options);


%compute t_hit
%obtain data from location
aux = 1e6*ones(length(options.x0),1);
inv = interval(-aux, aux);
tran_old = get(loc, 'transition');
reset = get(tran_old{1}, 'reset');
tran{1}=transition(h,reset,0,'a','b'); 
%create dummy location
loc = location('dummyLoc',0,inv,tran,contDynamics);  
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









function [tmin,tmax,t_total,RbeforeInt] = intersectionTime(obj,loc_old,h,options)

%initialize reachable set
contDynamics = get(loc_old,'contDynamics');

%positive direction
[tmin_pos,tmax_pos,t_total_pos,RbeforeInt_pos] = intersectionTime_oneDirection(h,obj,contDynamics,options);

%reverse dynamics
A = get(contDynamics,'A'); %get system matrix
B = get(contDynamics,'B'); %get input matrix
contDynamics_rev = linearSys('reverse_Sys',-A,-B); %initialize linear interval system;

%positive direction
[tmin_neg,tmax_neg,t_total_neg,RbeforeInt_neg] = intersectionTime_oneDirection(h,obj,contDynamics_rev,options);

%if negative direction is not empty 
if ~isempty(tmin_neg)
    if ~isempty(tmin_pos)
        tmin = -max(tmax_neg);
        tmax = max(tmax_pos);
        t_total = [];
        RbeforeInt = RbeforeInt_pos;
    else
        tmin = -tmax_neg;
        tmax = -tmin_neg;
        t_total = [];
        RbeforeInt = RbeforeInt_neg;
    end
else
    tmin = tmin_pos;
    tmax = tmax_pos;
    t_total = t_total_pos;
    RbeforeInt = RbeforeInt_pos;
end




% function [tmin,tmax,t_total,RbeforeInt] = intersectionTime_oneDirection(obj,contDynamics,h,options)
% 
% %refine
% options.timeStep = 0.1*options.timeStep;
% 
% %initialize reachable set
% [contDynamics,Rnext,options]=initReach(contDynamics,obj,options);
% 
% %is distance incearsing?
% delta_d = abs(h.c.'*center(Rnext.tp)-h.d) - abs(h.c.'*center(obj)-h.d);
% if delta_d > 0
%     dist_incr = 1; 
% else
%     dist_incr = 0; 
% end
% 
% %check if still in front of guard set
% intersect = zonoIntersect(h, Rnext.ti);
% 
% if ~dist_incr || intersect
% 
%     noIntersectionCounter = 0;
%     if intersect
%         RbeforeInt = obj;
%     else
%         %while reachable set still in front of guard set
%         while ~intersect
%             %set Rprev
%             Rprev = Rnext;
% 
%             %compute next reachable set
%             [Rnext, options]=reach(contDynamics, Rprev, options);
% 
%             %check if still in front of guard set
%             intersect = zonoIntersect(h, Rnext.ti);
% 
%             %increase counter
%             noIntersectionCounter = noIntersectionCounter + 1;
%         end
%         RbeforeInt = Rprev.tp;
%     end
% 
%     %obtain intersection duration
%     intersectionCounter = 0;
%     while intersect
%         %compute next reachable set
%         [Rnext, options]=reach(contDynamics, Rnext, options);
% 
%         %check if still in front of guard set
%         intersect = zonoIntersect(h, Rnext.ti);
% 
%         %increase counter
%         intersectionCounter = intersectionCounter + 1;
%     end
% 
%     %obtain t_min, t_max
%     tmin = 0;
%     tmax = intersectionCounter*options.timeStep;
%     t_total = (intersectionCounter + noIntersectionCounter - 1)*options.timeStep;
%     
% else
%     tmin = [];
%     tmax = [];
%     t_total = [];
%     RbeforeInt = [];
% end



function hIntermediate = intermediateHalfspaces(h_first, h_last, rotPoint, nrOfHalfspaces)


%obtain rotation matrix
%get dimension
dim = length(h_first.c);
%normalize normal vectors
n = h_first.c/norm(h_first.c);
newDir = h_last.c/norm(h_last.c);
%create mapping matrix
B(:,1) = n;
%find orthonormal basis for n, uVec
indVec = newDir - (newDir.'*n)*n;
B(:,2) = indVec/norm(indVec);
%complete mapping matrix B
if dim>2
    null(B(1:2,:)); %complete coding here
end
%compute angle between uVec and n
angle = acos(newDir.'*n);

%init
h = h_first;
delta_angle = angle;

for i = 1:nrOfHalfspaces
    
    delta_angle = delta_angle/2;
    
    %rotation matrix
    R = eye(dim);
    R(1,1) = cos(delta_angle);
    R(1,2) = -sin(delta_angle);
    R(2,1) = sin(delta_angle);
    R(2,2) = cos(delta_angle);
    %final rotation matrix
    rotMat = B*R*inv(B);

    %rotate halfspace
    h = rotMat*(h + (-rotPoint)) + rotPoint;

    %save halfspace
    hIntermediate{i} = h;
end
%rotate halfspace
h = rotMat*(h + (-rotPoint)) + rotPoint;

%save halfspace
hIntermediate{i+1} = h;


function hIntermediate = intermediateHalfspaces_linear(h_first, h_last, rotPoint, nrOfHalfspaces)


%obtain rotation matrix
%get dimension
dim = length(h_first.c);
%normalize normal vectors
n = h_first.c/norm(h_first.c);
newDir = h_last.c/norm(h_last.c);
%create mapping matrix
B(:,1) = n;
%find orthonormal basis for n, uVec
indVec = newDir - (newDir.'*n)*n;
B(:,2) = indVec/norm(indVec);
%complete mapping matrix B
if dim>2
    B(:,3:dim) = null(B(:,1:2).'); 
end
%compute angle between uVec and n
angle = acos(newDir.'*n);

%init
h = h_first;
delta_angle = angle/nrOfHalfspaces;

for i = 1:nrOfHalfspaces
    
    %rotation matrix
    R = eye(dim);
    R(1,1) = cos(delta_angle);
    R(1,2) = -sin(delta_angle);
    R(2,1) = sin(delta_angle);
    R(2,2) = cos(delta_angle);
    %final rotation matrix
    rotMat = B*R*inv(B);

    %rotate halfspace
    h = rotMat*(h + (-rotPoint)) + rotPoint;

    %save halfspace
    hIntermediate{i} = h;
end




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