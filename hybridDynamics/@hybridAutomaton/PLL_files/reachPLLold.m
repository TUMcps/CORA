function [obj,Rtotal,location] = reachPLLold(obj,options)
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
% Last revision: ---

%------------- BEGIN CODE --------------

%load data from options
tFinal=options.tFinal; %final time

%initialize variables
tInter=options.tStart; %intermediate time at transitions
nextLoc=options.startLoc; %next location
Rinter=options.R0; %intermediate reachable set at transitions
location=nextLoc; %location history vector
Rremain=[];

if options.startLoc==1
    pickOne = 1;
end

%while there exists a next reachable set
while ~isempty(Rinter)
    %specify new input
    options.U=options.Uloc{nextLoc};
    options.u=options.uLoc{nextLoc};
    options.uTrans=options.uLocTrans{nextLoc};

    %compute reachable set within a location
    options.t=tInter; %<-- for analysis purposes
    %set time step
    if nextLoc==1
        options.timeStep=options.timeStepFast;
    else
        options.timeStep=options.timeStepSlow;
    end
    
    %plot phase difference
    phaseDiff=interval([0 0 0 1 -1]*Rinter)
    
    [TP,R,loc,Rjump,Rcont]=reach(obj.location{nextLoc},tInter,Rinter,options);
    
    phaseDiff=interval([0 0 0 1 -1]*Rcont{end})
    
    %obtain vector of minimal times for reachin the guard sets
    tMin=TP.tMin;
    
    TP
    
    if ~isempty(tMin)
        %save results in a list
        if ~isempty(Rjump)
            for i=1:length(loc)
                %save variables of remaining paths
                ind=length(Rremain)+1;
                Rremain{ind}=Rjump{i};
                tRemain(ind)=tMin(i); 
                locRemain(ind)=loc(i);
            end
        end

        locRemain

        %search for jobs 
        anyOne=any(locRemain==1);
        allOne=all(locRemain==1);
        %picking state machine
        if pickOne && ~anyOne
            pickOne=0;
        end
        if ~pickOne && allOne
            pickOne=1;
        end    
        


        if tMin>4
            disp('stop')
        end
        
        if pickOne
            nextInd=find(locRemain==1,1);
        else
            %unify sets of locations 2,3
            for iLoc=2:3
                index = find(locRemain==iLoc);
                if ~isempty(index)
                    %compute enclosing set
                    if length(index)==2
                        %enclose parallel sets
                        try
                        Rencl = enclose(Rremain{index(1)}, Rremain{index(2)});
%                         phaseDiff=interval([0 0 0 1 -1]*Rremain{index(1)})
%                         phaseDiff=interval([0 0 0 1 -1]*Rremain{index(2)})
%                         phaseDiff=interval([0 0 0 1 -1]*Rencl)
                        catch
                            disp('merging error');
                        end
                    elseif length(index)==1
                        Rencl = Rremain{index};
                    end
                    %compute new minimum time
                    tMin=min(tRemain(index));
                    %delete seperate sets
                    Rremain(index)=[];
                    tRemain(index)=[]; 
                    locRemain(index)=[];
                    %add new sets
                    ind=length(Rremain)+1;
                    Rremain{ind}=Rencl;
                    tRemain(ind)=tMin; 
                    locRemain(ind)=iLoc;
                end
            end
            nextInd=find(locRemain~=1,1);
        end

    else
        %choose new index
        nextInd=1;
    end
    
    if ~isempty(Rremain)
        %determine variables of chosen path
        Rinter=Rremain{nextInd};
        tInter=tRemain(nextInd);
        nextLoc=locRemain(nextInd);

        %delete sets
        Rremain(nextInd)=[];
        tRemain(nextInd)=[]; 
        locRemain(nextInd)=[];
    else
        Rinter=[];
    end
    
    %store results
    index=length(location);
    location(end+1)=nextLoc;    
    Rtotal.OT{index}=R; 
    Rtotal.T{index}{1}=R{end}; %overapproximate time point polytope by last time interval polytope
    TPtotal{index}=TP;
    RcontTotal.OT{index}=Rcont;
    RcontTotal.T{index}{1}=Rcont{end};
end

%save results to object structure
obj.result.reachSet.R=Rtotal;
obj.result.reachSet.location=location;
obj.result.reachSet.TP=TPtotal;
obj.result.reachSet.Rcont=RcontTotal;


%------------- END OF CODE --------------