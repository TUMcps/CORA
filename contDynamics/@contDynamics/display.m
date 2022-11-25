function display(obj)
% display - Displays a continuous dynamics object
%
% Syntax:  
%    display(obj)
%
% Inputs:
%    obj - contDynamics object
%
% Outputs:
%    ---
%
% Example: 
%    cd=contDynamics('test function',[1 2],1,3);
%    display(cd);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 02-May-2007
% Last update: 17-October-2007
% Last revision: ---

%------------- BEGIN CODE --------------

%display name and id
disp(['Continuous dynamics "',obj.name,'"']);

% display number of states
disp(['number of states: ', num2str(obj.dim)]);

% display number of inputs
disp(['number of inputs: ' num2str(obj.nrOfInputs)]);

% display number of outputs
disp(['number of outputs: ' num2str(obj.nrOfOutputs)]);

%------------- END OF CODE --------------