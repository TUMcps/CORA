function [result] = in(obj1,obj2)
% in - determines if elements of a zonotope obj2 are in a halfspace
% obj1 in an overapproximate way
%
% Syntax:  
%    [result] = in(obj1,obj2)
%
% Inputs:
%    obj1 - halfspace object
%    obj2 - zonotope object
%
% Outputs:
%    result - 1/0 if zonotope is in, or not
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      14-June-2016
% Last update:  27-July-2016
% Last revision:---

%------------- BEGIN CODE --------------

%the zonotope penetrates the halfspace if part of the projection is
%negative

%1st case: obj2 is a single zonotope
if ~iscell(obj2)
    %projection
    projZ = obj1.c'*obj2 - obj1.d;
    projIH = interval(projZ);
    
    %minimum negative
    result = (infimum(projIH) <= 0);

%2nd case: obj2 is a cell array
else
    for i=1:length(obj2)
        %code to be entered...
    end
end
%------------- END OF CODE --------------