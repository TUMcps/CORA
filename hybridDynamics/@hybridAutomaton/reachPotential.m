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
% Last revision: ---

%------------- BEGIN CODE --------------

%load data from options
tFinal=options.tFinal; %final time

%initialize variables
tInter=options.tStart; %intermediate time at transitions
nextLoc=options.startLoc; %next location
%Rinter{1}=options.R0; %intermediate reachable set at transitions
Rinter=options.R0; %intermediate reachable set at transitions
openIndices=[]; %contains indices of unfinished reachable sets
location=nextLoc; %location history vector


%while there exists a next reachable set
while ~isempty(Rinter)
    %specify new input
    options.U=options.Uloc{nextLoc};
    options.u=options.uLoc{nextLoc};
    options.uTrans=options.uLocTrans{nextLoc};
    options.timeStep=options.timeStepLoc{nextLoc};

    %compute reachable set within a location
    [TP,R,loc,Rjump,Rcont]=reach(obj.location{nextLoc},tInter,Rinter,options);
    %obtain vector of minimal times for reachin the guard sets
    tMin=TP.tMin;
    loc=setdiff(loc,options.finalLoc); %remove final location
    
    if ~isempty(loc) && (tMin(1)<(tFinal*(1-1e-6))) 
        %determine variables of first path
        Rinter=Rjump{1};
        tInter=tMin(1);
        nextLoc=loc(1);
        
        %in case there are further paths
        for i=2:length(loc)
            %add open index
            if ~exist('tRemain')
                openIndices=1;
            else
                openIndices(end+1)=length(tRemain)+1;
            end
            %save variables of remaining paths
            ind=openIndices(end);
            Rremain{ind}=Rjump{i};
            tRemain(ind)=tMin(i); 
            locRemain(ind)=loc(i);
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
            %delete first open index
            openIndices(1)=[];
        else
            %no further open reachable sets
            Rinter=[];
        end
%         %save final set
%         for i=1:length(R)
%             if i==1
%                 Rfinal=R{i}{end};
%             else
%                 Rfinal=Rfinal & R{i}{end};
%             end
%         end
%         obj.result.final.R{finalIndex}=Rfinal;
%         obj.result.final.location(finalIndex)=nextLoc;
%         finalIndex = finalIndex+1;
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

% zon2bRed=RcontTotal.OT; %<--store to assess halfspace conversion methods 
% save zon2bRed zon2bRed;

%------------- END OF CODE --------------