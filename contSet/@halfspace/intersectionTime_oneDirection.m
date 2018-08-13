function [tmin,tmax,t_total,Rlast] = intersectionTime_oneDirection(obj,R0,contDynamics,options)
% intersectionTime_oneDirection - computes the time to intersect the
% halfspace in the given flow direction
%
% Syntax:  
%    [tmin,tmax,t_total,Rlast] = intersectionTime_oneDirection(obj,R0,contDynamics,options)
%
% Inputs:
%    obj - halfspace object
%    R0 - initial set
%    contDynamics - continuous dynamics
%    options - options struct
%
% Outputs:
%    tmin - minimum time
%    tmax - maximum time
%    t_total - time elapsed for complete intersection starting at R0
%    Rlast - last reachable set before first intersection with halfspace
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      05-September-2013
% Last update:  25-August-2016
% Last revision:---

%------------- BEGIN CODE --------------

% %test
% figure

%refine
options.timeStep = 0.1*options.timeStep;

%initialize reachable set
[Rnext,options]=initReach(contDynamics,R0,options);

%is distance incearsing?
delta_d = obj.c.'*(center(Rnext.tp) - center(R0));
if delta_d > 0
    dist_incr = 1; 
else
    dist_incr = 0; 
end

%check if still in front of guard set
intersect = zonoIntersect(obj, Rnext.ti);

if ~dist_incr || intersect

    noIntersectionCounter = 0;
    if intersect
        Rlast = R0;
    else
        %while reachable set still in front of guard set
        while ~intersect
            %set Rprev
            Rprev = Rnext;

            %compute next reachable set
            [Rnext, options]=post(contDynamics, Rprev, options);
            
%             %test
%             plot(Rnext.ti);

            %check if still in front of guard set
            intersect = zonoIntersect(obj, Rnext.ti);

            %increase counter
            noIntersectionCounter = noIntersectionCounter + 1;
        end
        Rlast = Rprev.tp;
    end

    %obtain intersection duration
    intersectionCounter = 0;
    while intersect
        %compute next reachable set
        [Rnext, options]=post(contDynamics, Rnext, options);
        
%         %test
%         plot(Rnext.ti,[1 2],'r');

        %check if still in front of guard set
        intersect = zonoIntersect(obj, Rnext.ti);

        %increase counter
        intersectionCounter = intersectionCounter + 1;
    end

    %obtain t_min, t_max
    tmin = 0;
    tmax = intersectionCounter*options.timeStep;
    t_total = (intersectionCounter + noIntersectionCounter - 1)*options.timeStep;
    
%     %test
%     plot(obj);
    
else
    tmin = [];
    tmax = [];
    t_total = [];
    Rlast = [];
end

%------------- END OF CODE --------------