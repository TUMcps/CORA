function [handle] = eventFcn(obj)
% eventFcn - returns the handle of the event function for a location
%
% Syntax:  
%    [handle] = eventFcn(obj)
%
% Inputs:
%    obj - location object
%
% Outputs:
%    handle - event function handle
%
% Example: 
%
% Other m-files required: eventFcn
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      07-May-2007 
% Last update:  06-June-2011
%               17-October-2013
% Last revision:---

%------------- BEGIN CODE --------------


%standard syntax for event functions as needed in Matlab simulations
function [value,isterminal,direction] = f(t,x)
    %get result of invariant events
    if ~isempty(obj.invariant)
        [value,isterminal,direction] = eventFcn(obj.invariant,x,1);
    else
        value = [];
        isterminal = [];
        direction = [];
    end
    %retrieve system dimension
    dim=length(value);
    %get result of guard events
        for i=1:length(obj.transition)
            %check if guard is a halfspace
             [resValue,resIsterminal,resDirection] = eventFcn(obj.transition{i},x);
             eventLength = length(resValue);
             indices = dim+1:dim+eventLength;
             value(indices,1) = resValue;
             isterminal(indices,1) = resIsterminal;
             direction(indices,1) = resDirection;       
             dim = dim+eventLength;
        end
end

%return function handle
handle = @f;

end

%------------- END OF CODE --------------
