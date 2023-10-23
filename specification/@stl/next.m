function res = next(obj,time)
% next - next-operator for Signal Temporal Logic
%
% Syntax:
%    res = next(obj,time)
%
% Inputs:
%    obj - logic formula (class stl)
%    time - scalar representing the time when the formula is active
%
% Outputs:
%    res - resulting stl formula (class stl)
%
% Example: 
%    x = stl('x',2);
%    eq = next(x(1) < 5,1.2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Authors:       Niklas Kochdumper, Benedikt Seidl
% Written:       09-November-2022 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % check input arguments
    if ~isa(obj,'stl') ||  ~obj.logic
        throw(CORAerror('CORA:notSupported',...
                      'This operation is not supported for stl objects!'));
    end
    
    % construct resulting stl object
    res = obj;
    
    res.type = 'next';
    res.lhs = obj;
    res.rhs = [];
    res.from = time;
    res.to = [];
    res.id = [];
    res.temporal = true;
end

% ------------------------------ END OF CODE ------------------------------
