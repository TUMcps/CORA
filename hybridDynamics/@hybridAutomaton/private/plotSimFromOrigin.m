function plotSimFromOrigin(obj,options)
% plot - plots simulation results of a hybrid automaton where the initial
% state is translated to the origin
%
% Syntax:  
%    plot(obj)
%
% Inputs:
%    obj - hybrid automaton object
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

% Author:       Matthias Althoff
% Written:      03-September-2009
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%load data from object structure
t=obj.result.simulation.t;
x=obj.result.simulation.x;
location=obj.result.simulation.location;

dim1=options.projectedDimensions(1);
dim2=options.projectedDimensions(2);

%plot trajectory segments
%         for i=1:length(t)
%             plot(t{i},x{i}(:,1),'-r');
%             hold on
%             plot(t{i},x{i}(:,2),'-g');
%         end
%         %plot switching times
%         for i=1:length(t)
%             vline(t{i}(end),':c',['t_{',num2str(i),'}']);
%         end
%         %plot switching planes
%         hline(1.9,'-k',['switching plane']);
%         hline(0.1,'-k',['switching plane']);
% 
% 
%         xlabel('time');
%         ylabel('x_1, x_2');
%         legend('x_1','x_2');    

%plot phase plane
for i=1:length(t)
    %plot(x{i}(:,dim1),x{i}(:,dim2),'-g');
    %plot(x{i}(:,dim1),x{i}(:,dim2),'-k');
    plot(x{i}(:,dim1)-x{i}(1,dim1),x{i}(:,dim2)-x{i}(1,dim2),'Color',0.6*[1 1 1]);
    hold on
end
xlabel(['x_',num2str(dim1)]);
ylabel(['x_',num2str(dim2)]);


%------------- END OF CODE --------------