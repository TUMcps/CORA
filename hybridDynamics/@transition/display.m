function display(obj)
% display - Displays a transition object
%
% Syntax:  
%    display(obj)
%
% Inputs:
%    obj - transition object
%
% Outputs:
%    ---
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 20-September-2007
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%display guard
display(obj.guard);

% display synchronization label
disp(strcat('Synchronization label:',obj.synchLabel));

%display target
disp(['target:',num2str(obj.target)]);

%------------- END OF CODE --------------