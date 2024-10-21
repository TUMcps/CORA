% TODO: rename as Minkowski difference
function int = minkDiff(obj,sub,varargin)
% minkDiff - compute the Minkowski difference of two STL intervals:
%            I1 - I2 = I <-> I + I2 \subseteq I1
%
% Syntax:
%    int = minkDiff(obj,sub)
%    int = minkDiff(obj,sub,type)
%
% Inputs:
%    obj - stlInterval object
%    sub - stlInterval object
%    type - type of computation ('exact' or 'inner')
%
% Outputs:
%    int - interval object after Minkowski difference
%
% Example: 
%    I1 = stlInterval(0,1);
%    I2 = stlInterval(0.5,1.5);
%
%    I = minkDiff(I1,I2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Florian Lercher
% Written:       06-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
type = setDefaultValues({'exact'},varargin);

% check input arguments
inputArgsCheck({{obj,'att','stlInterval'};
                {sub,'att',{'stlInterval'}};
                {type,'str',{'exact','inner'}}});

if isemptyobject(sub)
    int = stlInterval(0,inf,true,false);
    return;
end
if isemptyobject(obj)
    int = stlInterval();
    return;
end
lb = obj.lower - sub.lower;
ub = obj.upper - sub.upper;
lc = obj.leftClosed || ~sub.leftClosed;
rc = obj.rightClosed || ~sub.rightClosed;

% clamp lower bound to 0 if it is negative
if lb < 0
    lb = 0;
    lc = true;
end

int = stlInterval(lb,ub,lc,rc);

% ------------------------------ END OF CODE ------------------------------
