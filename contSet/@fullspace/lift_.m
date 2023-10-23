function fs = lift_(fs,N,dims)
% lift_ - lifts a full-dimensional space onto a higher-dimensional space
%    case R^0: undefined
%
% Syntax:
%    fs = lift_(fs,N,dims)
%
% Inputs:
%    fs - fullspace object
%    N - dimension of the higher-dimensional space
%    dims - states of the high-dimensional space that correspond to the
%          states of the low-dimensional space
%
% Outputs:
%    fs - lifted fullspace
%
% Example: 
%    fs = fullspace(4);
%    val = lift(fs,6,[1,2,5,6]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/lift

% Authors:       Mark Wetzlinger
% Written:       06-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
% for dims: 'size',[1,fs.dimension] ... but does not work always due to
% transposition of vector
inputArgsCheck({{fs,'att','fullspace'};
                {N,'att','numeric',{'scalar','nonnegative','integer','>=',fs.dimension}};
                {dims,'att','numeric',{'vector','nonnegative','integer'}}});

fs.dimension = N;

% ------------------------------ END OF CODE ------------------------------
