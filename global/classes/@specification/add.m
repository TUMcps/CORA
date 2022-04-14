function obj = add(obj1,obj2)
% add - joins two specification objects
%
% Syntax:  
%    obj = add(obj1,obj2)
%
% Inputs:
%    obj1 - specification object
%    obj2 - specification object
%
% Outputs:
%    obj - resulting specification object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: specification

% Author:       Niklas Kochdumper
% Written:      29-May-2020             
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    obj = repelem(obj1(1,1),size(obj1,1)+size(obj2,1),1);
    cnt = 1;
    
    for i = 1:size(obj1,1)
        obj(cnt,1) = obj1(i,1);
        cnt = cnt + 1;
    end
    
    for i = 1:size(obj2,1)
        obj(cnt,1) = obj2(i,1);
        cnt = cnt + 1;
    end

%------------- END OF CODE --------------