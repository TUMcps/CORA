function [list] = indexList(obj)
% indexList - returns a list that relates active events to guards; an
% invariant is defined as 0, other numbers refer to the guard number
%
% Syntax:  
%    [list] = indexList(obj)
%
% Inputs:
%    obj - location object
%
% Outputs:
%   list - list of guard numbers, whereas the list position refers to the
%   event
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 07-May-2007 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%get indices of the invariant
eq = get(obj.invariant,'equations');
list(1:eq) = 0;

%get indices of the guards
for i=1:length(obj.transition)
    eq = get(obj.transition{i},'equations');
    list((end+1):(end+eq)) = i;
end

%------------- END OF CODE --------------