function start_carReach(~)
% start_carReach - starts the function carReach
%
% Syntax:  
%    start_carReach()
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
[fileName,pathName] = uigetfile('','Load fArray');

% set model initialization
modelInitialization = @initCar;

% generate probabilistic model
probModel = carReach(fileName,pathName,modelInitialization);

% obtain interaction with other car
ThetaC = interaction(fileName,pathName,modelInitialization);

%save results to file
[file,path] = uiputfile('*.mat','Save probabilistic model as');
cd(path);
save(file,'probModel','ThetaC');


%------------- END OF CODE --------------