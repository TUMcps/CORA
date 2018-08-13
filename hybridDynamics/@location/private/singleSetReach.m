function [TP,R,activeGuards,Rguard,Rcont] = singleSetReach(obj,tStart,R0,options)
% singleSetReach - computes the reachable set of the system within a 
% location, detects the guard set that is hit and computes the reset - for
% a single initial set.
%
% Syntax:  
%    [TP,R,nextLoc,Rjump,Rcont] = singleSetReach(obj,tStart,R0,options)
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
%    nextLoc - next location
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

% Author: Matthias Althoff
% Written: 07-May-2007 
% Last update:  17-August-2007
%               26-March-2008
%               11-October-2008
%               21-April-2009
%               06-July-2009
%               10-August-2010
%               21-October-2010
%               25-October-2010
%               27-July-2016
%               30-July-2016
%               08-February-2017
%               20-April-2017 (option to intersect invariant made explicit)
% Last revision: ---

%------------- BEGIN CODE --------------

%load data from options
tFinal=options.tFinal;

%initialize set counter and time
iSet=1;
t=tStart;
inInv=1;

%if a trajectory should be tracked
if isfield(options,'uTransVec')
    options.uTrans = options.uTransVec(:,1);
end

%initialize reachable set
[Rnext,options]=initReach(obj.contDynamics,R0,options);

% if iscell(Rnext.tp)
%     Rnext.tp = Rnext.tp{1};
%     Rnext.ti = Rnext.ti{1};
% end

Rcont.OT{iSet}=Rnext.ti; 
Rcont.T{iSet}=Rnext.tp;

%increment time and set counter
t = t+options.timeStep;
iSet = iSet+1; 


%while reachable set still in invariant
while inInv && (t<tFinal)
    
    %if a trajectory should be tracked
    if isfield(options,'uTransVec')
        options.uTrans = options.uTransVec(:,iSet);
    end
    
    %compute reachable set within location
    [Rnext, options]=post(obj.contDynamics, Rnext, options);
    
%     %work-around for cell-structured reachable sets (it is assumed that no splitting occured)
%     if iscell(Rnext.tp)
%         Rnext.tp = Rnext.tp{1};
%         Rnext.ti = Rnext.ti{1};
%     end
    
    Rcont.OT{iSet}=Rnext.ti;     
    Rcont.T{iSet}=Rnext.tp;
    %increment time and set counter
    t = t+options.timeStep;
    iSet = iSet+1; 

    if iscell(Rnext.ti)
        for i=1:length(Rnext.ti)
            inInvVec(i) = in(obj.invariant,Rnext.ti{i}); %check if reachable set is in invariant
        end
        inInv = any(inInvVec);
    else
        inInv=in(obj.invariant,Rnext.ti); %check if reachable set is in invariant 
    end
end

% %work-around for cell-structured reachable sets (it is assumed that no splitting occured)
% for i = 1:length(Rcont.OT)
%     if iscell(Rcont.OT{i})
%         Rcont.OT{i} = Rcont.OT{i}{1};
%     end
% end

%check if reachable sets potentially intersect with guards
[guards,setIndices]=potInt(obj,Rcont.OT);

%compute intersections of overapproximated reachable sets
[Rguard,activeGuards,minIndex,maxIndex] = guardIntersect(obj,guards,setIndices,Rcont.OT,options);

%remove reachable parts outside the invariant
if isfield(options,'intersectInvariant') && options.intersectInvariant==1
    R = potOut(obj,Rcont.OT,minIndex,maxIndex,options);
else
    R=Rcont;
end
        

%if guard set has been hit
if ~isempty(minIndex)  
    %obtain time point struct TP
    deltaMin=(minIndex-1)*options.timeStep;
    deltaMax=(maxIndex-1)*options.timeStep;
    
    TP.tStart=tStart;
    TP.tMin=tStart+deltaMin;
    TP.tMax=tStart+deltaMax;

%if no guard set has been hit    
else
    TP.tStart=tStart;
    TP.tMin=t;
    TP.tMax=t;

    Rguard{1}=[];
    activeGuards=[];
end
%------------- END OF CODE --------------