function display(intMat)
% display - Displays an intervalMatrix object on the command window
%
% Syntax:  
%    display(intMat)
%
% Inputs:
%    intMat - intervalMatrix object
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

% Author:       Matthias Althoff
% Written:      18-June-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%display dimension, generators
disp('dimension: ');
disp(intMat.dim);
%display left and right limits
disp('left limit: ');
disp(infimum(intMat.int));
disp('right limit: ');
disp(supremum(intMat.int));

%------------- END OF CODE --------------