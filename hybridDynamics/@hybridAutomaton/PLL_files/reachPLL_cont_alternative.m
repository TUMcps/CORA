function [R] = reachPLL_cont_alternative(obj,options)
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

%initialize
totalCycles = 0;
t_cycle = 1/options.sys.c(5);
dim = length(options.sys.A);
Rgoal = options.Rgoal;
Rinit = options.Rinit;
cycles = 10;

activeRegions = 1:length(Rgoal);

for i=1:40
    
    %obtain linear map
    totalCycles = totalCycles + cycles;
    t = totalCycles*t_cycle;
    linMap = expm(options.sys.A*t) + 0.01*eye(dim);
    linMap_inv = pinv(linMap);
    
    %update options.W
    options.W{1} = eye(dim);
    options.W{2} = linMap;

    %compute reachable set starting at each goal region
    for iRegion = activeRegions
        if iRegion==1
            options.overflow=1;
        else
            options.overflow=0;
        end
        %compute fixed number of cycles
        Rfinal{iRegion} = reachPLL_cycles(Rinit{iRegion}, cycles, options);

        phaseDiff=interval([0 0 0 1 -1]*Rfinal{iRegion})

    end
    
    
    %intersect reachable sets with goal regions
    for iRegion = activeRegions
        for iSet = activeRegions
            %intersect
            Rintersect{iRegion}{iSet} = Rgoal{iRegion} & Rfinal{iSet};
            Rintersect_map{iRegion}{iSet} = linMap_inv*Rintersect{iRegion}{iSet};
        end
    end
    
    %unify
    for iRegion = activeRegions
        IH = interval(Rintersect{iRegion}{activeRegions(1)});
        IH_map = interval(Rintersect_map{iRegion}{activeRegions(1)});
        for iSet = activeRegions(2:end)
            %unify intersected reachable sets
            try
            IH = IH | interval(Rintersect{iRegion}{iSet});
            catch
                a=1
            end
            IH_map = IH_map | interval(Rintersect_map{iRegion}{iSet});
        end
        %convert to zonotope
        %if ~isempty(get(IH,'intervals')) && ~isempty(get(IH_map,'intervals'))
        try
            Zinit{1}=zonotope(IH);
            Zinit{2}=linMap*zonotope(IH_map);
            Rinit{iRegion} = zonotopeBundle(Zinit);
        %else
        catch
            %delete region
            ind = find(activeRegions==iRegion);
            activeRegions(ind) = [];
        end
        
%         %plot
%         dims=[1 4];
%         for iSet = activeRegions(1:end)
%             plot(Rintersect{iRegion}{iSet},dims);
%         end
%         plot(Rinit{iRegion},dims,'r');
    end

    

%     %update variables
%     totalCycles = totalCycles + cycles;
%     t = totalCycles*t_cycle;
%     linMap = expm(options.sys.A*t) + 0.01*eye(dim);
%     W{1} = eye(dim);
%     W{2} = linMap;
%     options.W = W;
% 
%     %intersect reachable sets with goal regions
%     for iRegion = 1:length(Rgoal)
%         for iSet = 1:length(Rgoal)
%             %intersect
%             Rtmp = Rgoal{iRegion} & Rfinal{iSet};
%             P{iRegion}{iSet} = polytope(Rtmp,options);
%         end
%     end
% 
%     for iRegion = 1:length(Rgoal)
%         %unify intersected reachable sets
%         Rinit{iRegion}=enclosePolytopes(P{iRegion},options);
%     end
    

Rstore{i}=Rinit;
end
disp('stop');

function [Zbundle] = enclosePolytopes(P,options)

%obtain new vertices
V=[];
sets=length(P);

%obtain vertices of each polytope
for iSet=1:sets
    %check if polytope is empty
    if ~isempty(P{iSet})
        Vpartial=vertices(P{iSet});
        V=[V,get(Vpartial,'V')];
    end
end
V=vertices(V);

%direction matrices
W=options.W;

for i=1:length(W)

    %exponential matrix transformation
    Z{i} = W{i}*zonotope(interval(pinv(W{i})*V));
end

Zbundle=zonotopeBundle(Z);

% figure;
% hold on
% dims=[1 2];
% plot(Zbundle,dims);
% for iStep=1:length(P)
%     plot(P{iStep},dims,'g');
% end
% 
% figure;
% hold on
% dims=[2 3];
% plot(Zbundle,dims);
% for iStep=1:length(P)
%     plot(P{iStep},dims,'g');
% end
% 
% figure;
% hold on
% dims=[4 5];
% plot(Zbundle,dims);
% for iStep=1:length(P)
%     plot(P{iStep},dims,'g');
% end


%------------- END OF CODE --------------