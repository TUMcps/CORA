function res = test_capsule_plot
% test_capsule_plot - unit test function of plot; this function aims
%    to go through many variations of input arguments
%    note: only run-time errors checked, go through manually to check for bugs
%
% Syntax:  
%    res = test_capsule_plot
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

% Author:       Mark Wetzlinger
% Written:      04-August-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% instantiate capsule
C = capsule([1;1;2],[0;1;-0.5],0.5);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plot(C);
    
    % two arguments: object, dimensions
    plot(C,1);
    plot(C,[1,2]);
    plot(C,[2,3]);
    
    % three arguments: object, dimensions, linespec
    plot(C,[1,2],'r+');
    
    % three arguments: object, dimensions, NVpairs
    plot(C,[1,2],'LineWidth',2);
    plot(C,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(C,[1,2],'FaceColor','r','LineWidth',2);
    plot(C,[1,2],'FaceColor','r','LineWidth',2,'EdgeColor',[.6 .6 .6]);
    
    % close figure
    close;
catch
    close;
    res = false;
end

%------------- END OF CODE --------------
