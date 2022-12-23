function R = add(R1,R2,varargin)
% add - joins two reachSet objects
%
% Syntax:  
%    R = add(R1,R2)
%    R = add(R1,R2,parent)
%
% Inputs:
%    R1 - reachSet object
%    R2 - reachSet object
%    parent - index of the parent for the root of the reachSet object obj2
%
% Outputs:
%    R - resulting reachSet object
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

% one of the objects is empty
if isempty(R1)
    R = R2;
    return
elseif isempty(R2)
    R = R1;
    return
end

% parse input arguments
parent = setDefaultValues({0},varargin);

% check input arguments
inputArgsCheck({{R1,'att',{'reachSet'},{''}};
                {R2,'att',{'reachSet'},{''}};
                {parent,'att',{'numeric'},{'integer','nonnegative','scalar'}}});

% add objects together
R = repelem(R1(1,1),size(R1,1)+size(R2,1),1);
cnt = 1;

for i = 1:size(R1,1)
    R(cnt,1) = R1(i,1);
    cnt = cnt + 1;
end

for i = 1:size(R2,1)
    R(cnt,1) = R2(i,1);
    if R(cnt,1).parent == 0
        R(cnt,1).parent = parent;
    else
        R(cnt,1).parent = R(cnt,1).parent + size(R1,1);
    end
    cnt = cnt + 1;
end

%------------- END OF CODE --------------