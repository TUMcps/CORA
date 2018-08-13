function [TP,R,nextLoc,blockedLoc,Rjump,Rcont] = reach(obj,tStart,R0,blockedLoc,options)
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
%               08-August-2011
%               11-September-2013
%               31-July-2016
%               19-August-2016
% Last revision: ---

%------------- BEGIN CODE --------------


%set up TPtmp
TP.tMin = [];
TP.tMax = [];

%set up variables
nextLoc = [];
Rjump = [];
Rguard_noInt = [];

if options.isHyperplaneMap
    % get results from hyperplane mapping
    %[TP,R,activeGuards,Rguard,Rguard_noInt,Rcont] = hyperplaneMap(obj,tStart,R0,blockedLoc,options);
    [TP,R,activeGuards,Rguard,Rguard_noInt,Rcont] = hyperplaneMap_noInv(obj,tStart,R0,blockedLoc,options);
    Rguard = Rguard(activeGuards);
    Rguard_noInt = Rguard_noInt(activeGuards);
    TP.tMin = TP.tMin(activeGuards);
else

    %get results from reachability analysis
    [TP,R,activeGuards,Rguard,Rcont] = singleSetReach(obj,tStart,R0,options);

end

%compute jumps, next locations
counter = 1;
TPsaved = TP;
blockedLoc = [];

for i=1:length(activeGuards)

    %obtain guard number
    iGuard = activeGuards(i);

    %does no intersection set exist?
    if ~isempty(Rguard_noInt)
        for iSet = 1:length(Rguard_noInt{i})
            if ~isempty(Rguard_noInt{i}{iSet})
                Rjump{counter} = reset(obj.transition{iGuard},Rguard_noInt{i}{iSet});  
                nextLoc{counter,1} = obj.id;
                blockedLoc{counter,1} = get(obj.transition{iGuard},'target');
                TP.tMin(counter) = TPsaved.tMin(i);

                %update counter
                counter = counter + 1;
            end
        end
    end

    %compute reset and new locations of active guards
    for iSet = 1:length(Rguard{i})
        Rjump{counter} = reset(obj.transition{iGuard},Rguard{i}{iSet});  
        nextLoc{counter,1}=get(obj.transition{iGuard},'target');
        blockedLoc{counter,1} = 0;
        TP.tMin(counter) = TPsaved.tMin(i);

        %update counter
        counter = counter + 1;
    end
end


%------------- END OF CODE --------------