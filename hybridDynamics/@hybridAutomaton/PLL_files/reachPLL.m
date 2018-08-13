function [obj,Rtotal,location] = reachPLL(obj,options)
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
%               30-July-2016
% Last revision: ---

%------------- BEGIN CODE --------------

%load data from options
tFinal=options.tFinal; %final time

%initialize variables
tInter=options.tStart; %intermediate time at transitions
nextLoc=options.startLoc; %next location
Rinter=options.R0; %intermediate reachable set at transitions
location=nextLoc; %location history vector
Rstore=[];

cycleCount=0;
lastIncrease=0;
iscontained=[];

%mapping matrix
M=[-1 0 0 0 0;...
   0 -1 0 0 0;...
   0 0 -1 0 0;...
   0 0 0 0 1;...
   0 0 0 1 0];


%while there exists a next reachable set
while ~isempty(Rinter) && ~(interval(Rinter)<=options.Rgoal) || (cycleCount<=30)

    %specify new input
    options.U=options.Uloc{nextLoc};
    options.u=options.uLoc{nextLoc};
    options.uTrans=options.uLocTrans{nextLoc};

    %compute reachable set within a location
    options.t=tInter; %<-- for analysis purposes
    %set time step
    if nextLoc==1
        options.timeStep=options.timeStepFast;
        %count cycle
        cycleCount=cycleCount+1
    else
        options.timeStep=options.timeStepSlow;
    end
    
    %plot phase difference
    phaseDiff=interval([0 0 0 1 -1]*Rinter)
    vi=interval([1 0 0 0 0]*Rinter)
    vp1=interval([0 1 0 0 0]*Rinter)
    
    [TP,R,loc,Rjump,Rcont]=reach(obj.location{nextLoc},tInter,Rinter,options);
    
    phaseDiff=interval([0 0 0 1 -1]*Rcont{end})
    
    %obtain vector of minimal times for reachin the guard sets
    tMin=TP.tMin;
    
    TP
    
    if ~isempty(loc)
        
        %store result of location that is not considered
        removeInd=find(loc==options.unconsideredLocation,1);
        otherInd=find(loc~=options.unconsideredLocation,1);
        if ~isempty(removeInd)
            lastIncrease
%             if isempty(Rstore)
%                 Rstore=interval(Rjump{removeInd});
%             else
%                 IHtmp = interval(Rjump{removeInd});
%                 if ~(IHtmp<=Rstore)
%                     Rstore=Rstore | IHtmp;
%                     lastIncrease=cycleCount;
%                 end
%             end
            
            if isempty(Rstore)
                Rstore=Rjump{removeInd};
                Pstore=polytope(Rstore);
            else
                Pjump=polytope(Rjump{removeInd});
                if ~(Pstore.iscontained(Pjump))
                    Ptmp{1}=Pstore;
                    Ptmp{2}=Pjump;
                    options.W{2}=eye(Ptmp{1}.dim);
                    Rstore=enclosePolytopes(obj.location{loc(removeInd)},Ptmp,options);
                    %slightly increase Rstore
                    Rstore=enlarge(Rstore,[1.01, 1.01, 1.01, 1.001, 1.01]);
                    Pstore=polytope(Rstore);
                    lastIncrease=cycleCount;
                end
            end            
            
            
%             %enlarge Rjump
%             %Rjump{otherInd}=enlarge(Rjump{otherInd},[1.02, 1.02, 1.02, 1.001, 1.02]);
%             Rjump{otherInd}=enlarge(Rjump{otherInd},[1.2, 1.2, 1.2, 1.02, 1.2]);
%             %check for containment
%             Plarge = polytope(Rjump{otherInd});
%             Psmall = polytope(M*Rjump{removeInd});
%             
%             if Plarge.iscontained(Psmall)
%                 iscontained(end+1)=1;
%             else
%                 iscontained(end+1)=0;
%             end
            %remove unconsidered location
            loc(removeInd)=[];
            Rjump(removeInd)=[];
            tMin(removeInd)=[];
        end

        disp('locations:');
        loc
        nextLoc=loc(1);
        %determine variables of chosen path
        Rinter=Rjump{1};
        tInter=tMin(1);
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
obj.result.reachSet.Rstore=Rstore;
obj.result.reachSet.Rjump=Rjump;
obj.result.reachSet.cycleCount=cycleCount;
%obj.result.reachSet.iscontained=iscontained;
obj.result.reachSet.lastIncrease=lastIncrease;



%------------- END OF CODE --------------