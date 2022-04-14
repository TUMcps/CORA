function [id,msg] = errDimMismatch()
% errDimMismatch - standardized error message format for dimension mismatch
%    of objects in set operations
%
% Syntax:  
%    errStruct = errDimMismatch()
%
% Inputs:
%    ---
%
% Outputs:
%    id - error identifier
%    msg - error message
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: noops

% Author:        Mark Wetzlinger
% Written:       12-March-2021
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

id = 'CORA:dimensionMismatch';
msg = 'Dimension mismatch in input arguments!';

%------------- END OF CODE --------------