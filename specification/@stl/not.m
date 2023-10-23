function res = not(obj)
% not - overloads the ~ operator representing negation for stl objects
%
% Syntax:
%    res = not(obj)
%
% Inputs:
%    obj - logic formula (class stl)
%
% Outputs:
%    res - resulting stl formula (class stl)
%
% Example: 
%    x = stl('x',2);
%    eq = ~(x(1) < 5)
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
    
    res.type = '~';
    res.lhs = obj;
    res.rhs = [];
    res.from = [];
    res.to = [];
    res.id = [];
end

% ------------------------------ END OF CODE ------------------------------
