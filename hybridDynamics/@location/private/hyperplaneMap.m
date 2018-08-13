function [TP,R,activeGuards,Rguard,Rguard_noInt,Rcont] = hyperplaneMap(obj,tStart,R0,blockedLoc,options)
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
% Written:      08-August-2011
% Last update:  10-September-2013
% Last revision:---

%------------- BEGIN CODE --------------

%load data from options
tFinal=options.tFinal;

%initialize set counter and time
iSet=1;
t=tStart;


%initialize reachable set
[obj.contDynamics,Rnext,options]=initReach(obj.contDynamics,R0,options);
Rcont.OT{iSet}=Rnext.ti; 
Rcont.T{iSet}=Rnext.tp; 
Rfull.T{iSet}=Rnext.tp;
%increment time and set counter
t = t+options.timeStep;
iSet = iSet+1; 

%check if still in invariant
inInv = in(obj.invariant,Rnext.ti);

%while reachable set still COMPLETELY in invariant
while inInv && (t<tFinal)
    
    %compute next reachable set
    [Rnext, options]=reach(obj.contDynamics, Rnext, options);

    Rcont.OT{iSet}=reduce(Rnext.ti,'girard',options.zonotopeOrder);
    Rcont.T{iSet}=reduce(Rnext.tp,'girard',options.zonotopeOrder);
    Rfull.T{iSet}=Rnext.tp;
    %increment time and set counter
    t = t+options.timeStep;
    iSet = iSet+1;
    
    if mod(iSet,5000)==0
        iSet
        length(get(Rnext.ti,'Z'))
    end
    
    if options.isHybrid
        %check if still in invariant
        inInv = in(obj.invariant,Rnext.ti);
    end
end

%check if system is hybrid
if options.isHybrid

    %detect intersections of overapproximated reachable sets
    [activeGuards,minIndex,maxIndex] = hyperplaneGuardIntersect(obj,Rcont.OT); 
    
    %remove reachable parts outside the invariant
    %[R] = potOut(obj,Rcont.OT,options);
    R=Rcont.OT;  
    
    %remove blocked locations
    iGuard = 1;
    while iGuard <= length(activeGuards)
        
        %obtain guard nr
        guardNr = activeGuards(iGuard);
        
        %obtain next location
        nextLoc = get(obj.transition{guardNr},'target');
        
        %remove blocked guards
        if ismember(nextLoc, blockedLoc)
            activeGuards(iGuard) = [];
        else
            iGuard = iGuard + 1;
        end
    end
    
    %work around for HSCC 14b
    if length(activeGuards)>1
        activeGuards(2:end) = [];
    end

    
    %compute Rguard
    for iGuard = 1:length(activeGuards)

%         disp('GUARD INTERSECT TIME');
%         tic
        
        %obtain guard Nr
        guardNr = activeGuards(iGuard);

        tmin = 0;
        tmax = (maxIndex(guardNr)+1 - minIndex(guardNr))*options.timeStep;
        A = get(get(obj,'contDynamics'),'A');
        B = get(get(obj,'contDynamics'),'B');
        guard = get(obj.transition{activeGuards(iGuard)},'guard');
        
        % set before intersection
        firstInd = minIndex(guardNr)-1;
        if firstInd == 0
            RbeforeInt = R0;
        else
            RbeforeInt = Rfull.T{firstInd};
        end
        
        % set after intersection
        lastInd = maxIndex(guardNr)+1;
        if lastInd-1 == length(Rfull.T)
            RafterInt = Rfull.T{lastInd-1};
        else
            RafterInt = Rfull.T{lastInd};
        end

        %reduce zonotope
        RbeforeInt=reduce(RbeforeInt,'girard',options.zonotopeOrder);
        RafterInt=reduce(RafterInt,'girard',options.zonotopeOrder);

        %compute guard projection
        partialIntersection = ~zonoIn(guard, RafterInt);
        
        if partialIntersection
            options.maxProjectionError = 0.05;
            %options.maxProjectionError = 0.06;
            %options.maxProjectionError = 0.1;
            %options.maxProjectionError = 0.5;
        else
            options.maxProjectionError = 0.05;
            %options.maxProjectionError = 0.1;
            %options.maxProjectionError = 0.3;
            %options.maxProjectionError = 0.5;
        end
        
        %set counter
        counter = 1;
        counter_noInt = 1;
        R_new{1} = RbeforeInt;
        
        while counter <= length(R_new)
            %set split flag
            splitFlag = 1;
            
            while splitFlag
                [R_guard, R_guard_noInt, splitFlag, RbeforeInt] = guardProjection(R_new{counter},guard,A,B,obj,partialIntersection,options);

                if splitFlag
                    %split
%                     %normalize normal vectors
%                     n = guard.c/norm(guard.c);
%                     newDir = flowVec/norm(flowVec);
%                     %find orthonormal vector
%                     orthDir = newDir - (newDir.'*n)*n;
%                     orthDir = orthDir/norm(orthDir);
%                     
%                     %projected flow
%                     projectedFlow = (flowVec.'*n)*n + (flowVec.'*orthDir)*orthDir;
                    
                    %R_aux = split(RbeforeInt,guard.c,'bundle');
                    R_aux = split(RbeforeInt,guard.c);

                    %delete current R_new
                    R_new(counter) = [];
                    %determine number of sets
                    nrOfSets = length(R_new);
                    %insert split sets
                    R_new{nrOfSets + 1} = R_aux{1};
                    R_new{nrOfSets + 2} = R_aux{2};
                else
                    %store guard intersections
                    for iSet = 1:length(R_guard)
                        R_new{counter} = R_guard{iSet};
                        counter = counter + 1;
                        
                        %plot
                        if options.debug
                            dim = length(center(R_new{counter - 1}));
                            Rred = reduce(R_new{counter - 1}, 'girard', 10);
                            Rred = Rred + zonotope([zeros(dim,1), 1e-4*eye(dim)]);
                            plot(Rred, [1 2], 'g');
                        end
                    end
                    %store set with no intersection
                    R_new_noInt{counter_noInt} = R_guard_noInt;
                    counter_noInt = counter_noInt + 1;
                end
            end
        end
        
%         %unify
%         R_union2 = R_new{1};
%         for iSet = 2:length(R_new)
%             R_union2 = enclose(R_union2,R_new{iSet});
%         end
%         if length(R_new) > 1
%             R_union = or(R_new{1}, R_new(2:end));
%         else
%             R_union = R_new{1};
%         end

        %add uncertain input
        optionsTmp = options;
        optionsTmp.timeStep = tmax;
        Rinput = errorSolution(obj.contDynamics,B*options.U,optionsTmp);

%         toc
        %set that has not intersected
        for i = 1:length(R_new_noInt)
            if ~isempty(R_new_noInt{i})
                R_new_noInt{i} = R_new_noInt{i} + Rinput;
                %reduce zonotope
                Rguard_noInt{guardNr}{i}=reduce(R_new_noInt{i},'girard',options.zonotopeOrder);
            else
                Rguard_noInt{guardNr}{i}=[];
            end
        end

%         R_union = R_union + Rinput;
%         %reduce zonotope
%         Rguard{guardNr}{1}=reduce(R_union,'girard',options.zonotopeOrder);
        
        for i = 1:length(R_new)
            R_new{i} = R_new{i} + Rinput;
            %reduce zonotope
            Rguard{guardNr}{i}=reduce(R_new{i},'girard',options.zonotopeOrder);
        end

        
    end
else
    minIndex=[];
    R=Rcont;
end

%if guard set has been hit
if ~isempty(activeGuards)  
    %obtain time point struct TP
    deltaMin = minIndex*options.timeStep;
    deltaMax = maxIndex*options.timeStep;
    
    TP.tStart = tStart;
    TP.tMin = tStart + deltaMin;
    TP.tMax = tStart + deltaMax;

%if no guard set has been hit    
else
    TP.tStart = tStart;
    TP.tMin = t;
    TP.tMax = t;

    Rguard{1} = [];
    Rguard_noInt{1} = [];
    activeGuards = [];
end
%------------- END OF CODE --------------