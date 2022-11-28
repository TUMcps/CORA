function display(obj)
% display - Displays a contDynamics object on the command window
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
%    sys = contDynamics('test function',[1 2],1,3);
%    display(sys);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      02-May-2007
% Last update:  17-October-2007
%               19-June-2022
%               23-November-2022 (TL: dispInput)
% Last revision:---

%------------- BEGIN CODE --------------


% disp input if necessary
dispInput(inputname(1))

%display name and id
disp("Continuous dynamics: '" + obj.name + "'");

% display number of states
disp("  number of states: " + obj.dim);

% display number of inputs
disp("  number of inputs: " + obj.nrOfInputs);

% display number of outputs
disp("  number of outputs: " + obj.nrOfOutputs);

fprintf(newline);

%------------- END OF CODE --------------