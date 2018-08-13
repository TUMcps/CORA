function [R] = deleteRedundantSets(R,Rold,options)
% deleteRedundantSets - delete reachable sets that are already covered by
% other sets
%
% Syntax:  
%    [R] = deleteRedundantSets(R,Rold,options)
%
% Inputs:
%    R - reachable sets
%    Rold - reachable sets of previous time steps
%    options - options for the computation of the reachable set
%
% Outputs:
%    R - reachable sets
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      29-June-2009
% Last update:  03-February-2011
% Last revision:---

%------------- BEGIN CODE --------------


%set reduction method
%redMethod='methC';
redMethod='PP';

%increase internal count
try
    R.internalCount=Rold.internalCount+1;
catch
    %for first run
    R.internalCount=3;
end

%if internal count='some constant'
if R.internalCount==options.reductionInterval
    %reset internal count
    R.internalCount=1;
    %overapproximate reachable set of time point(!) by parallelpipeds and
    %save them as polytopes
    R.P=[];
    for i=1:length(R.tp)
        R.tp{i}=reduce(R.tp{i},redMethod);
        %generate polytope
        Ptmp=polytope(R.tp{i});
        R.P{i}=mptPolytope(Ptmp.C,Ptmp.d);
    end    
elseif R.internalCount==2
    %intersect each reachable set with each previous reachable set
    for iNewSet=1:length(R.tp)
        %approximate new set of time points by parallelpiped 
        R.tp{iNewSet}=reduce(R.tp{iNewSet},redMethod);
        %generate mpt polytope
        Ptmp=polytope(R.tp{iNewSet});
        Pnew{iNewSet}=mptPolytope(Ptmp.C,Ptmp.d);
    end
    %initialize Pcut
    Pcut=Pnew;
    %intersection with previous time step
    for iNewSet=1:length(Pcut)
        %intersect with previous sets
        for iOldSet=1:length(Rold.P)           
            Pcut{iNewSet}=Pcut{iNewSet}\Rold.P{iOldSet};
        end
    end
    
    %reset iChecked counter
    iChecked=1;
    
    %figure;
    %hold on
    
    %intersection with actual time step
    for iNewSet=1:length(Pcut)
        %intersect with actual sets
        for iOtherSet=1:length(Pnew)  
            if iOtherSet~=iNewSet
                Pcut{iNewSet}=Pcut{iNewSet}\Pnew{iOtherSet};
            end
        end
        %is polytope empty?
        if ~isempty(Pcut{iNewSet})
            Rnew{iChecked}=R.tp{iNewSet};
            iChecked=iChecked+1;
            %plot(Rnew{iChecked-1});
        else
            disp('canceled!!');
        end    
    end
    
    %copy only ckecked reachable sets
    R.tp=[];
    R.tp=Rnew;
end


%------------- END OF CODE --------------