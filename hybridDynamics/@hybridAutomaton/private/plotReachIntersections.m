function plotReachIntersections(obj,options)
% plotReachIntersections - plots intersections of reachable sets of a 
% hybrid automaton
%
% Syntax:  
%    plotReachIntersections(obj,options)
%
% Inputs:
%    obj - hybrid automaton object
%    options - options struct
%
% Outputs:
%    none
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 11-October-2008
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

hold on

%load data from object structure
R=obj.result.reachSet.R.OT;
loc=obj.result.reachSet.location;
dim=options.projectedDimensions;
nrOfParallelSets=[];


%plot each reachable set
for i=1:length(R)
    %determine number of time steps
    timeSteps=inf;
    for iSet=1:length(R{i})
        steps=length(R{i}{iSet});
        if steps<timeSteps
            timeSteps=steps;
        end
    end    
    
    for iStep=1:timeSteps
        Pint=[];
        for iIntersect=1:length(R{i})
            %if reachable set  is a zonotope
            if strcmp('zonotope',class(R{i}{iIntersect}{iStep}))
                Rred=project(R{i}{iIntersect}{iStep},dim);
                %Rred=reduce(Rred,'girard',6);
                Rred=reduce(Rred,'girard',1.5);
                Pproj=polytope(Rred);
            %if reachable set is a polytope
            elseif strcmp('polytope',class(R{i}{iIntersect}{iStep})) %polytope
                %project polytope
                Pproj=projection(R{i}{iIntersect}{iStep},dim); 
            end
            
            if isempty(Pint)
                Pint=Pproj;
            else
                %Pint=Pint&Pproj;
                if iIntersect==1
                    Pint=Pproj;
                end
            end
        end
        %plot(Pint);
        V=extreme(Pint)';
        V=vertices(V);
        %V=extreme(Pproj)';
        
        
        plot(V,options.plotType);
    end
end


%------------- END OF CODE --------------