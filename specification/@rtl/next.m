function res = next(obj,time)
% next - next-operator for Reachset Temporal Logic
%
% Syntax:
%    res = next(obj,time)
%
% Inputs:
%    obj - logic formula (class rtl)
%    time - scalar representing the time when the formula is active
%
% Outputs:
%    res - resulting rtl formula (class rtl)
%
% Example: 
%    x = stl('x',2)
%    eq = next(rtl(x(1) < 5),1.0)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Authors:       Niklas Kochdumper
% Written:       09-November-2022 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{obj,'att','rtl'}});

% ensure that time is a multiple of 0.5
if ~withinTol(mod(time,0.5),0)
    throw(CORAerror('CORA:wrongValue','second',...
    'Argument "time" for next operator has to be integer or multiple of 0.5!'));
end

% construct resulting rtl object
res = obj;

res.type = 'next';
res.lhs = obj;
res.rhs = [];
res.time = time;

% ------------------------------ END OF CODE ------------------------------
