function start_carReach_reach(~)
% start_carReach_reach - starts the function carReach_reach
%
% Syntax:  
%    start_carReach_reach(~)
%
% Inputs:
%    no
%
% Outputs:
%    no (result is saved as file)
%
% Example: 
%
% 
% Author:       Matthias Althoff
% Written:      31-July-2017
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------


%load fArray to determine segment length of road 
[fileName,pathName] = uigetfile();

% set model initialization
modelInitialization = @initCar;

% generate probabilistic model
probModel = carReach_reach(fileName,pathName,modelInitialization);

%save results to file
[file,path] = uiputfile('*.mat','Save probabilistic model as');
cd(path);
save(file,'probModel');


%------------- END OF CODE --------------