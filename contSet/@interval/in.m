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

TOL = 1e-12;
res = false;

% point in interval containment
if isnumeric(obj2)
    
    % account for numerical jitter
    obj2lb = obj2 + TOL;
    obj2ub = obj2 - TOL;
    
    if all([all(obj1.inf <= obj2lb), all(obj1.sup >= obj2ub)])
        res = true;
    end

% interval in interval containment
elseif isa(obj2,'interval')
    
    % check for dimension mismatch
    if dim(obj1) ~= dim(obj2)
        [id,msg] = errDimMismatch();
        error(id,msg);
    end

    if all(obj1.sup >= obj2.sup) && all(obj1.inf <= obj2.inf)
        res = true;
    end

% other set in interval containment
else

    res = in(mptPolytope(obj1),obj2);

end

%------------- END OF CODE --------------