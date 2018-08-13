function [TP,R,guardNr,Rguard,Rguard_noInt,Rcont] = hyperplaneMap_noInv(obj,tStart,R0,blockedLoc,options)
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
%               19-August-2016
% Last revision:---

%------------- BEGIN CODE --------------

%load data from options
tFinal=options.tFinal;

%initialize set counter and time
iSet=1;
t=tStart;
guardNr = [];

%initialize reachable set
[Rnext,options]=initReach(obj.contDynamics,R0,options);
Rcont.OT{iSet}=Rnext.ti; 
Rcont.T{iSet}=Rnext.tp; 
Rfull.T{iSet}=Rnext.tp;


%check if no guard is hit
noGuard = isempty(hyperplaneIntersectCheck(obj,Rnext.ti));

%while reachable set still COMPLETELY in invariant
while noGuard && (t<tFinal)
    
    %compute next reachable set
    [Rnext, options]=post(obj.contDynamics, Rnext, options);
    
    %increment time and set counter
    t = t+options.timeStep;
    iSet = iSet+1; 

    Rcont.OT{iSet}=reduce(Rnext.ti,'girard',options.zonotopeOrder);
    Rcont.T{iSet}=reduce(Rnext.tp,'girard',options.zonotopeOrder);
    Rfull.T{iSet}=Rnext.tp;

    %check if still in invariant
    noGuard = isempty(hyperplaneIntersectCheck(obj,Rnext.ti));
end

%remove reachable parts outside the invariant
%[R] = potOut(obj,Rcont.OT,options);
R=Rcont.OT;

%check if system is hybrid
if ~noGuard  
    
    %obtain guard nr
    guardNr = hyperplaneIntersectCheck(obj,Rnext.ti);
    currGuard = get(obj.transition{guardNr},'guard');

    %compute Rguard
    firstInd = iSet - 1;
    
    %Rfirst
    RbeforeInt = Rcont.T{firstInd};
    
    %continue until fully in guard
    maxGuardIntersection = 2;
    while (length(hyperplaneIntersectCheck(obj, Rnext.ti)) < maxGuardIntersection) && ...
            ~zonoIn(currGuard, Rnext.ti) && (t<tFinal)
    
        %compute next reachable set
        [Rnext, options]=post(obj.contDynamics, Rnext, options);

        Rcont.OT{iSet}=reduce(Rnext.ti,'girard',options.zonotopeOrder);
        Rcont.T{iSet}=reduce(Rnext.tp,'girard',options.zonotopeOrder);
        Rfull.T{iSet}=Rnext.tp;
        
        %plot(Rcont.OT{iSet});
        
        %increment time and set counter
        t = t+options.timeStep;
        iSet = iSet+1;    
    end
    
    % set after intersection
    lastInd = iSet - 1;
    RafterInt = Rfull.T{lastInd};

    %reduce zonotope
    RbeforeInt=reduce(RbeforeInt,'girard',options.zonotopeOrder);
    RafterInt=reduce(RafterInt,'girard',options.zonotopeOrder);

    %compute guard projection
    partialIntersection = ~zonoIn(currGuard, RafterInt);
    
    %obtain data for guard intersection
    A = get(get(obj,'contDynamics'),'A');
    B = get(get(obj,'contDynamics'),'B');

    %set counter
    counter = 1;
    counter_noInt = 1;
    R_new{1} = RbeforeInt;
        
    while counter <= length(R_new) && ~isempty(R_new{1})
        %set split flag
        splitFlag = 1;

        while splitFlag
            if isa(currGuard,'constrainedHyperplane')
                currGuard = currGuard.h;
            end
            [R_guard, R_guard_noInt, splitFlag, RbeforeInt] = guardProjection(R_new{counter},currGuard,A,B,obj,partialIntersection,options);

            if splitFlag

                %R_aux = split(RbeforeInt,guard.c,'bundle');
                R_aux = split(RbeforeInt,currGuard.c);

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
                    if ~isempty(R_guard{iSet})
                        R_new{counter} = R_guard{iSet};

                        %plot
                        if isfield(options,'debug') && options.debug
                            dim = length(center(R_new{counter - 1}));
                            Rred = reduce(R_new{counter - 1}, 'girard', 10);
                            Rred = Rred + zonotope([zeros(dim,1), 1e-4*eye(dim)]);
                            plot(Rred, [1 2], 'g');
                        end
                    end
                    counter = counter + 1;
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

    %TODO: move computation to error E;add uncertain input
    optionsTmp = options;
    tmax = options.timeStep*(lastInd - firstInd);
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

%if guard set has been hit
if ~noGuard
    %obtain time point struct TP
    deltaMin = firstInd*options.timeStep;
    deltaMax = lastInd*options.timeStep;
    
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
end
%------------- END OF CODE --------------
