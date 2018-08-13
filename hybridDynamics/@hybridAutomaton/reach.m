function [obj,Rtotal,location] = reach(obj,options)
% reach - computes the reachable set of a hybrid automaton
%
% Syntax:  
%    [obj] = reach(obj,options)
%
% Inputs:
%    obj - hybrid automaton object
%    options - options for simulation, reachability analysis of systems
%
% Outputs:
%    obj - hybrid automaton object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      07-May-2007 
% Last update:  16-August-2007
%               26-March-2008
%               07-October-2008
%               21-April-2009
%               20-August-2013
%               30-October-2015
%               22-August-2016
% Last revision: ---

%------------- BEGIN CODE --------------

%load data from options
tFinal = options.tFinal; %final time

%initialize variables
tInter = options.tStart; %intermediate time at transitions
nextLoc = options.startLoc; %next location
nextBlockedLoc = [];
Rinter = options.R0; %intermediate reachable set at transitions
openIndices = []; %contains indices of unfinished reachable sets
tRemain = []; %contains indices of unfinished reachable sets
location = nextLoc; %location history vector


%while there exists a next reachable set
while ~isempty(Rinter)
    
    %specify new input
    options.U=options.Uloc{nextLoc};
    options.u=options.uLoc{nextLoc};
    if isfield(options,'uLocTransVec')
        options.uTransVec = options.uLocTransVec{nextLoc};
    else
        options.uTrans=options.uLocTrans{nextLoc};
    end
    options.timeStep=options.timeStepLoc{nextLoc};
    
    %obtain factors for initial state and input solution
    for i=1:(options.taylorTerms+1)
        r = options.timeStepLoc{nextLoc};
        options.factor(i)= r^(i)/factorial(i);    
    end
    
    if ~(nextLoc==options.finalLoc)
        % convert the set representation of the guard sets if necessary
        obj.location{nextLoc} = convGuard(obj.location{nextLoc},options);
        
        %compute reachable set within a location
        [TP,R,loc,blockedLoc,Rjump,Rcont]=reach(obj.location{nextLoc},tInter,Rinter,nextBlockedLoc,options);
        
        loc = cell2mat(loc);
        blockedLoc = cell2mat(blockedLoc);
        
        %obtain vector of minimal times for reachin the guard sets
        tMin=TP.tMin;
    else
        loc = [];
    end
        

    if ~isempty(loc) && (tMin(1)<tFinal) 
        %determine variables of first path
        Rinter=Rjump{1};
        tInter=tMin(1);
        nextLoc=loc(1);
        nextBlockedLoc=blockedLoc(1);
        
        %in case there are further paths
        for i=2:length(loc)
            %add open index
            openIndices(end+1)=length(tRemain)+1;

            %save variables of remaining paths
            ind=openIndices(end);
            Rremain{ind}=Rjump{i};
            tRemain(ind)=tMin(i); 
            locRemain(ind)=loc(i);
            locBlockedRemain(ind)=blockedLoc(i);
        end

    %if previous path is finished
    else
        %if there are still open indices
        if ~isempty(openIndices)
            %load variables of unsolved path
            ind=openIndices(1);
            Rinter=Rremain{ind};
            tInter=tRemain(ind);
            nextLoc=locRemain(ind);
            nextBlockedLoc=locBlockedRemain(ind);
            %delete first open index
            openIndices(1)=[];
        else
            %no further open reachable sets
            Rinter=[];
        end
    end 
    
    %store results
    index=length(location);
    location(end+1)=nextLoc;  
    try
        Rtotal.OT{index}=R.OT;
        Rtotal.T{index}{1}=R.T{end}; %overapproximate time point polytope by last time interval polytope
    catch
        Rtotal.OT{index}=R; 
        Rtotal.T{index}{1}=R{end}; 
    end
    TPtotal{index}=TP;
    RcontTotal.OT{index}=Rcont.OT;
    RcontTotal.T{index}=Rcont.T;
    Rencl{index}=Rjump;
    
end

%save results to object structure
obj.result.reachSet.R=Rtotal;
obj.result.reachSet.location=location(1:end-1);
obj.result.reachSet.TP=TPtotal;
obj.result.reachSet.Rcont=RcontTotal;
obj.result.reachSet.Rencl=Rencl;

%------------- END OF CODE --------------