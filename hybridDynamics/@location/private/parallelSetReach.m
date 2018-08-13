function [TP,R,Rcont,activeGuards,Rguard] = parallelSetReach(obj,tStart,R0,options)
% reach - computes the reachable set of the system within a location, 
% detects the guard set that is hit and computes the reset
%
% Syntax:  
%    [TP,R,nextLoc,Rjump,Rcont] = reach(obj,tStart,R0,options)
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
%    Rcont - reachable set due to continuous evolution
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
% Written:      07-May-2007 
% Last update:  17-August-2007
%               26-March-2008
%               21-April-2009
%               06-July-2009
%               24-July-2009
%               31-July-2009
%               11-August-2010
%               22-October-2010
%               25-October-2010
% Last revision: ---

%------------- BEGIN CODE --------------


    
parallelSets=length(R0);
TPpartial=cell(1,parallelSets);
R=cell(1,parallelSets);
activeGuardsPartial=cell(1,parallelSets);
RunorderedGuard=cell(1,parallelSets);
Rcont=cell(1,parallelSets);

%parfor (iSet=1:parallelSets)
for iSet=1:parallelSets
    %get results from reachability analysis
    [TPpartial{iSet},R{iSet},activeGuardsPartial{iSet},RunorderedGuard{iSet},Rcont{iSet}] = singleSetReach(obj,tStart,R0{iSet},options);  
end

%update nextLoc
for iSet=1:length(R0)
    %check for parallel sets
    if length(activeGuardsPartial{iSet})>1
        disp('more than one next location...');
    end          
    %intersect next locations
    if iSet==1
        activeGuards=activeGuardsPartial{iSet};
    else
        activeGuards=intersect(activeGuards,activeGuardsPartial{iSet});
    end
end

activeGuards

%collect TP data    
%init
if ~isempty(activeGuards)
    %obtain final TP data
    TP.tStart=TPpartial{1}.tStart;
    %set up TPtmp
    TP.tMin(1:length(activeGuards))=tStart;
    TP.tMax(1:length(activeGuards))=inf;

    %collect TP data
    for iSet=1:length(TPpartial)
        %update TP
        for iGuard=1:length(activeGuards) 
            %obtain guard nr
            guardNr=activeGuards(iGuard);

            %update TP
            %tMin
            if TPpartial{iSet}.tMin(guardNr)>TP.tMin(iGuard)
                TP.tMin(iGuard)=TPpartial{iSet}.tMin(guardNr);
            end
            %tMax
            if TPpartial{iSet}.tMax(guardNr)<TP.tMax(iGuard)
                TP.tMax(iGuard)=TPpartial{iSet}.tMax(guardNr);
            end        
        end
    end

    TP

    %obtain min and max index
    minInd=round((TP.tMin-TP.tStart)/options.timeStep)+1;
    maxInd=round((TP.tMax-TP.tStart)/options.timeStep)+1;

    %next location index loop
    iGuard=1;
    while iGuard<=length(activeGuards)

        %obtain guard nr
        guardNr=activeGuards(iGuard);

        %init counter
        counter=1;
        iStep=minInd(iGuard);
        continueFlag=1;
        while (iStep<=maxInd(iGuard)) && continueFlag
        %for iStep=minInd:maxInd
            %initialize R jump
            Rguard{iGuard}{counter}=RunorderedGuard{1}{guardNr}{iStep};

            %init iSet
            iSet=2;
            while (iSet<=length(RunorderedGuard)) && continueFlag   
                %check for emptiness
                if isempty(RunorderedGuard{iSet}{guardNr}{iStep})
                    continueFlag=0;
                end
                %intersect polytopes
                Rguard{iGuard}{counter}=Rguard{iGuard}{counter} & RunorderedGuard{iSet}{guardNr}{iStep};
                %increase iSet
                iSet=iSet+1;
            end

            %if Rguard is empty, delete it
            if (iSet>length(RunorderedGuard)) && ~isempty(Rguard{iGuard}{counter})
                counter=counter+1;
            end
            %increase iStep
            iStep=iStep+1;
        end
        try
        %in case the next location is empty
        if isempty(Rguard{iGuard}{1})
            activeGuards(iGuard)=[];
        else
            iGuard=iGuard+1;
        end
        catch
            disp('error');
        end
    end
end
        




%------------- END OF CODE --------------