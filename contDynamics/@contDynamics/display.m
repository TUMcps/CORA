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

%display state IDs
disp(['state IDs: ', mat2str(obj.stateIDs)]); 

%display input IDs
disp(['input IDs: ', mat2str(obj.inputIDs)]); 

%display output IDs
disp(['output IDs: ', mat2str(obj.outputIDs)]); 

%display system dimension
disp(['state dimension: ', num2str(obj.dim)]);


%------------- END OF CODE --------------