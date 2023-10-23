function display(intMat)
% display - displays an intervalMatrix object on the command window
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
%    intMat = intervalMatrix([1 2 3; 2 3 1],[1 0 2; 0 1 1])
%    intMat = intervalMatrix()
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       18-June-2010
% Last update:   03-April-2023 (MW, add empty case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isempty(intMat)
    
    dispEmptyObj(intMat,inputname(1));

else

    % display dimension
    disp('dimension: ');
    disp(dim(intMat));
    
    % display lower and upper bounds
    disp('lower bound: ');
    disp(infimum(intMat.int));
    disp('upper bound: ');
    disp(supremum(intMat.int));

end

% ------------------------------ END OF CODE ------------------------------
