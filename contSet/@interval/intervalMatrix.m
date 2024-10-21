function intMat = intervalMatrix(I)
% intervalMatrix - determines if an interval intersects a set
%
% Syntax:
%    intMat = intervalMatrix(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    intMat - intervalMatrix object
%
% Example: 
%    I = interval([2 1; -3 4],[3 2; 1 5]);
%    intervalMatrix(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: intervalMatrix

% Authors:       Tobias Ladner
% Written:       10-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

intMat = intervalMatrix(center(I),rad(I));
    
end

% ------------------------------ END OF CODE ------------------------------
