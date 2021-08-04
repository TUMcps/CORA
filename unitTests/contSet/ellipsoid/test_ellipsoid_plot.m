function res = test_ellipsoid_plot
% test_ellipsoid_plot - unit test function of plot
%
% Syntax:  
%    res = test_ellipsoid_plot
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      27-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
load cases.mat E_c
for i=1:length(E_c)
    E1 = E_c{i}.E1;
    Ed1 = E_c{i}.Ed1;
    E0 = E_c{i}.E0;
    
    res = tryPlot(E1) && tryPlot(Ed1) && tryPlot(E0);
    
    if ~res
        break;
    end
    
end


if res
    disp([mfilename,' successful']);
else
    disp([mfilename,' failed']);
end

end

%-- helper
function res = tryPlot(E)
res = true;
try
    % try all variations in plotting
    h = figure;
    
    % one argument: object
    plot(E);
    
    % three arguments: object, dimensions, linespec
    plot(E,[1,2],'r+');
    
    % three arguments: object, dimensions, NVpairs
    plot(E,[1,2],'LineWidth',2);
    plot(E,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    
    % three arguments: object, dimensions, NVpair 'Filled'
    plot(E,[1,2],'Filled',true);
    plot(E,[1,2],'Filled',true,'LineWidth',2);
    if isFullDim(project(E,[1,2]))
        plot(E,[1,2],'Filled',true,'EdgeColor','k','FaceColor',[.8 .8 .8]);
    end
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(E,[1,2],'r','Filled',true,'LineWidth',2);
    if isFullDim(project(E,[1,2]))
        plot(E,[1,2],'r','Filled',true,'LineWidth',2,'EdgeColor',[.6 .6 .6]);
    end
    
    % close figure
    close(h);
catch
    close;
    res = false;
end
end
%------------- END OF CODE --------------
