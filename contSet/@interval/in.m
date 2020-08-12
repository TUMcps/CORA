function res = in(obj1,obj2)
% in - determines if obj2 is contained in obj1
%
% Syntax:  
%    res = in(obj1,obj2)
%
% Inputs:
%    obj1 - interval object
%    obj2 - contSet object or point
%
% Outputs:
%    res - 1/0 if obj2 is contained in obj1, or not
%
% Example: 
%    int1 = interval([-1;-2],[2;3]);
%    int2 = interval([0;0],[1;1]);
%
%    in(int1,int2)
%    in(int1,[1;2])
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/in

% Author:       Niklas Kochdumper
% Written:      16-May-2018
% Last update:  02-Sep-2019
% Last revision:---

%------------- BEGIN CODE --------------

res = 0;

% point in interval containment
if isnumeric(obj2)

    if all(obj1.sup >= obj2) && all(obj1.inf <= obj2)
        res = 1;
    end

% interval in interval containment
elseif isa(obj2,'interval')

    if all(obj1.sup >= obj2.sup) && all(obj1.inf <= obj2.inf)
        res = 1;
    end

% other set in interval containment
else

    res = in(mptPolytope(obj1),obj2);

end

%------------- END OF CODE --------------