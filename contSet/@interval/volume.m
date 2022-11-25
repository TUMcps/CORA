function V = volume(obj)
% volume - Computes volume of an interval
%
% Syntax:  
%    V = volume(obj)
%
% Inputs:
%    obj - interval object
%
% Outputs:
%    V - volume
%
% Example: 
%    I = interval([1; -1], [2; 1]);
%    V = volume(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      24-July-2016 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

if ~isempty(obj)
    V = prod(2*rad(obj));
else
    V = 0;
end

%------------- END OF CODE --------------