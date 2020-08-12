function obj = add(obj1,obj2)
% add - joins two simResult objects
%
% Syntax:  
%    obj = add(obj1,obj2)
%
% Inputs:
%    obj1 - simResult object
%    obj2 - simResult object
%
% Outputs:
%    obj - resulting simResult object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachSet

% Author:       Niklas Kochdumper
% Written:      29-May-2020             
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    if isempty(obj1)
         obj = obj2;
    elseif isempty(obj2)
         obj = obj1;
    else

        if isempty(obj1.loc) ~= isempty(obj2.loc)
           error('Objects are not compatible!'); 
        end

        obj = obj1;
        obj.x = [obj1.x; obj2.x];
        obj.t = [obj1.t; obj2.t];
        obj.loc = [obj1.loc; obj2.loc];
    end

%------------- END OF CODE --------------