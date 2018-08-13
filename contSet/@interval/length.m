function L = length(obj)
% length - Overloads the opertor that returns the length of the longest 
% array dimension
%
% Syntax:  
%    L = length(obj)
%
% Inputs:
%    obj - interval object 
%
% Outputs:
%    L - length of the largest array dimension in X.
%
% Example: 
%    a=interval([-1 1], [1 2]);
%    length(a)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      18-November-2015 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%returns length of infimum
L = length(obj.inf);

%------------- END OF CODE --------------