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
%    parent - (optional) index of the parent for the root of R2
%
% Outputs:
%    R - resulting reachSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachSet

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       29-May-2020             
% Last update:   06-June-2023 (MW, simplify)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% one of the objects is empty 
if (isscalar(R1) || isnumeric(R1)) && isemptyobject(R1)
    R = R2;
    return
elseif (isscalar(R2) || isnumeric(R2)) && isemptyobject(R2)
    R = R1;
    return
end

% parse input arguments
parent = setDefaultValues({0},varargin);

% check input arguments
inputArgsCheck({{R1,'att','reachSet'};
                {R2,'att','reachSet'};
                {parent,'att',{'numeric'},{'integer','nonnegative','scalar'}}});

% length of R1
R1_length = length(R1);

% indices where R2 parent is 0 -> overwritten by provided value for parent
idxZero = [R2.parent] == 0;

% loop of all non-zero parents in R2, increment by length of R1
for i=1:length(R2)
    if idxZero(i)
        R2(i).parent = parent;
    else
        R2(i).parent = R2(i).parent + R1_length;
    end
end

% concatenate objects
R = [R1;R2];

% ------------------------------ END OF CODE ------------------------------
