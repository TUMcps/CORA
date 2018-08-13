function res = test_requiredToolboxes
% test_requiredToolboxes - checks if required toolboxes are installed
%
% Syntax:  
%    res = test_requiredToolboxes
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% 
% Author:       Matthias Althoff
% Written:      15-September-2016
% Last update:  04-May-2018
% Last revision:---


%------------- BEGIN CODE --------------

%check if symbolic toolbox is available
res_partial(1) = license('test','Symbolic_Toolbox');
if (res_partial(1)==0)
    disp('(symbolic toolbox missing)');
end

%check if CORA and MPT are in the MATLAB path
p = path;
res_partial(2) = ~isempty(strfind(p,'contDynamics')); %one of the folders in CORA
res_partial(3) = ~isempty(strfind(p,'contSet')); %one of the folders in CORA
res_partial(4) = ~isempty(strfind(p,'hybridDynamics')); %one of the folders in CORA
res_partial(5) = ~isempty(strfind(p,'mpt')); %should be found if MPT toolbox is installed

res = all(res_partial);

%------------- END OF CODE --------------
