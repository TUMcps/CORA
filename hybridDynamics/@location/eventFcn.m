function han = eventFcn(loc)
% eventFcn - returns the handle of the event function for a location
%
% Syntax:
%    han = eventFcn(loc)
%
% Inputs:
%    loc - location object
%
% Outputs:
%    han - event function handle
%
% Example:
%    ---
%
% Other m-files required: eventFcn
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       07-May-2007 
% Last update:   06-June-2011
%                17-October-2013
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    %standard syntax for event functions as needed in Matlab simulations
    function [value,isterminal,direction] = f(t,x)
        %get result of invariant events
        if ~representsa_(loc.invariant,'emptySet',eps)
            [value,isterminal,direction] = eventFcn(loc.invariant,x,1);
        else
            value = [];
            isterminal = [];
            direction = [];
        end
        %retrieve system dimension
        n=length(value);
        %get result of guard events
        for i=1:length(loc.transition)
            %check if guard is a halfspace
            if ~isemptyobject(loc.transition(i))
                [resValue,resIsterminal,resDirection] = eventFcn(loc.transition(i),x);
                eventLength = length(resValue);
                indices = n+1:n+eventLength;
                value(indices,1) = resValue;
                isterminal(indices,1) = resIsterminal;
                direction(indices,1) = resDirection;       
                n = n+eventLength;
            end
        end
    end
    
    % return function handle
    han = @f;

end

% ------------------------------ END OF CODE ------------------------------
