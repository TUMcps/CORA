function res = nonlinearOperator(obj,operator)
% nonlinearOperator - apply nonlinear funtion such as sin() to stl object
%
% Syntax:
%    res = nonlinearOperator(obj,operator)
%
% Inputs:
%    obj - logic formula (class stl)
%    operator - string representing the operator
%
% Outputs:
%    res - resulting stl formula (class stl)
%
% Example: 
%    x = stl('x',2);
%    eq = nonlinearOperator(x(1),'sin')
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

    % check if input is correct
    if isa(obj,'stl') && ~obj.logic

        res = obj;

        res.lhs = obj;
        res.rhs = [];
        res.type = operator;
        
    else
        throw(CORAerror('CORA:notSupported',...
                      'This operation is not supported for stl objects!'));
    end
end

% ------------------------------ END OF CODE ------------------------------
