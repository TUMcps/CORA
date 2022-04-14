function obj = add(obj1,obj2,varargin)
% add - joins two reachSet objects
%
% Syntax:  
%    obj = add(obj1,obj2)
%    obj = add(obj1,obj2,parent)
%
% Inputs:
%    obj1 - reachSet object
%    obj2 - reachSet object
%    parent - index of the parent for the root of the reachSet object obj2
%
% Outputs:
%    obj - resulting reachSet object
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

    % parse input arguments
    parent = 0;
    if nargin >= 3
       parent = varargin{1}; 
    end

    % add objects together
    if isempty(obj1)
        obj = obj2;
    elseif isempty(obj2)
        obj = obj1;
    else
        
        obj = repelem(obj1(1,1),size(obj1,1)+size(obj2,1),1);
        cnt = 1;

        for i = 1:size(obj1,1)
            obj(cnt,1) = obj1(i,1);
            cnt = cnt + 1;
        end

        for i = 1:size(obj2,1)
            obj(cnt,1) = obj2(i,1);
            if obj(cnt,1).parent == 0
                obj(cnt,1).parent = parent;
            else
                obj(cnt,1).parent = obj(cnt,1).parent + size(obj1,1);
            end
            cnt = cnt + 1;
        end
    end

%------------- END OF CODE --------------