function res = finally(obj,time)
% finally - finally-operator for Signal Temporal Logic
%
% Syntax:  
%    res = finally(obj,time)
%
% Inputs:
%    obj - logic formula (class stl)
%    time - time interval (class interval)
%
% Outputs:
%    res - resulting stl formula (class stl)
%
% Example: 
%    x = stl('x',2)
%    eq = finally(x(1) < 5,interval(0.1,0.2))
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Author:       Niklas Kochdumper
% Written:      9-November-2022 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % check input arguments
    if ~isa(obj,'stl') || ~obj.logic
        throw(CORAerror('CORA:notSupported',...
                      'This operation is not supported for stl objects!'));
    end
    
    if ~isa(time,'interval') || ~all(size(time) == [1,1])
        throw(CORAerror('CORA:wrongValue',...
                        'Wrong format for input argument "time"!'));
    end

    % construct resulting stl object
    res = obj;
    
    res.type = 'finally';
    res.lhs = obj;
    res.rhs = [];
    res.temporal = true;
    res.from = infimum(time);
    res.to = supremum(time);
end

%------------- END OF CODE --------------