function plotCarInteraction(obj)
% plot - plots scenario, where two cars drive on the same lane; the
% positions, velocities and their distance is plotted
%
% Syntax:  
%    plotCarInteraction(obj,options)
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
% Written: 22-October-2007 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------


%load data from object structure
t=obj.result.simulation.t;
x=obj.result.simulation.x;
loc=obj.result.simulation.location;

figure;
%plot each reachable set
for i=1:length(t)
    if length(t{i})>1
        %positions
        subplot(3,1,1);
        plot(t{i},x{i}(:,1),'-r');
        hold on
        plot(t{i},x{i}(:,3),'-b');

        %velocities
        subplot(3,1,2);
        plot(t{i},x{i}(:,2),'-r');
        hold on
        plot(t{i},x{i}(:,4),'-b');    

        %distance (from bumper to bumper, including car length)
        subplot(3,1,3);
        plot(t{i},x{i}(:,5));
    end
end


%------------- END OF CODE --------------