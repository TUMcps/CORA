function res = example_commonocean2cora_ZAM
% example_commonocean2cora_ZAM - example for the conversion of a
%    CommonOcean XML-file scenario to CORA model
%
% Syntax:
%    res = example_commonocean2cora_ZAM
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

% Authors:       Bruno Maione
% Written:       16-February-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% convert model
[statObs,dynObs,x0,goalSet,waterways,shallows] = ...
    commonocean2cora('ZAM_Tutorial-2_1_T-1',false);

% visualization
figure; hold on;

% waterways
for i = 1:length(waterways)
    plot(waterways{i}.set, 'FaceColor',[.7 .7 .7]);
end

% shallows
for i = 1:length(shallows)
    plot(shallows{i}.set,'FaceColor','cyan');
end

% goal set and start point
plot(goalSet{1}.set,[1,2],'FaceColor','r','FaceAlpha',0.5);
plot(x0.x,x0.y,'.g','MarkerSize',20);

% dynamic obstacles
for i = 1:length(dynObs)
    plot(dynObs{i}.set,[1,2],'FaceColor','b','FaceAlpha',0.5);
end

% static obstacles
for i = 1:length(statObs)
    plot(statObs{i}.set, 'FaceColor', 'y','EdgeColor','none');
end

% set axis limits
xlim([-600,600]);
ylim([-600,600]);
axis equal

% successful
res = true;

% ------------------------------ END OF CODE ------------------------------
