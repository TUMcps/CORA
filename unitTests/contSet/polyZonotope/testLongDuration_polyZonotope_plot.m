function res = testLongDuration_polyZonotope_plot
% testLongDuration_polyZonotope_plot - unit test function of plot;
%    this function aims to go through many variations of input arguments
%    note: only run-time errors checked, manual bug check necessary
%
% Syntax:  
%    res = testLongDuration_polyZonotope_plot
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

% create random polyZonotope
c = rand(4,1)-0.5*ones(4,1);
G = rand(4,6)-0.5*ones(4,6);
ind = datasample(1:6,4,'Replace',false);
G(:,ind) = G(:,ind)./10;
Grest = rand(4,2)-0.5*ones(4,2);
expMat = [eye(4), round(rand(4,2)*5)];
pZ = polyZonotope(c,G,Grest,expMat);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plot(pZ);
    
    % two arguments: object, dimensions
    plot(pZ,[1,2]);
    plot(pZ,[2,3]);
    
    % three arguments: object, dimensions, linespec
    plot(pZ,[1,2],'r+');
    
    % three arguments: object, dimensions, NVpairs
    plot(pZ,[1,2],'LineWidth',2);
    plot(pZ,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    
    % three arguments: object, dimensions, NVpair 'Filled'
    plot(pZ,[1,2],'Filled',true);
    plot(pZ,[1,2],'Filled',true,'LineWidth',2);
    plot(pZ,[1,2],'Filled',true,'EdgeColor','k','FaceColor',[.8 .8 .8]);
    
    % three arguments: object, dimensions, NVpair 'Splits'
    plot(pZ,[1,2],'Splits',6);
    plot(pZ,[1,2],'Splits',6,'LineWidth',2);
    plot(pZ,[1,2],'Splits',6,'Filled',true,'EdgeColor','k','FaceColor',[.8 .8 .8]);
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(pZ,[1,2],'r','Filled',true,'LineWidth',2);
    plot(pZ,[1,2],'r','Filled',true,'LineWidth',2,'EdgeColor',[.6 .6 .6]);
    
    % close figure
    close;
catch
    close;
    res = false;
end

if res
    disp('test_plot successful');
else
    disp('test_plot failed');
end

%------------- END OF CODE --------------

