classdef stl
% stl - class representing a Signal Temporal Logic (STL) formula
%
% Syntax:
%    obj = stl(name)
%    obj = stl(name,num)
%    obj = stl(true)
%    obj = stl(false)
%
% Inputs:
%    name - name of the STL-variable
%    num - dimension of the variable
%
% Outputs:
%    obj - generated stl variable
%
% Example:
%    x = stl('x',2);
%    eq = until(x(1) <= 5,x(2) > 3 & x(1) <= 2, interval(0.1,0.2))
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: rtl

% Authors:       Niklas Kochdumper, Benedikt Seidl
% Written:       09-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    type = 'variable'       % type of the opeartor ('until', '&', etc.)
    lhs = []                % object forming left hand side of operator
    rhs = []                % object forming right hand side of operator
    from = []               % start of the time interval
    to = []                 % end of the time interval
    variables = {}          % list of variables present in the formula
    var = []                % name of the current variable
    temporal = false        % temporal (true) or non-temporal (false)
    logic = false           % propositional formula (true) or not (false)
    id = []                 % unique identifier for the formula
end
    
methods
    
    % class constructor
    function obj = stl(name,varargin)
        
        % catch the case where the input is a logic value
        if islogical(name)
            if name
                obj.type = 'true';
            else
                obj.type = 'false';
            end
            obj.variables = [];
            obj.logic = true;
            return;
        end

        % parse and check input arguments
        if ~ischar(name)
           throw(CORAerror('CORA:wrongInputInConstructor',...
                        'Input "name" has to be a string!'));
        end
        
        num = [];
        
        if nargin > 1
            
            num = varargin{1};
        
            if ~isa(num, 'atomicProposition') && (~isnumeric(num) || ...
                    ~isscalar(num) || floor(num) ~= ceil(num))
                throw(CORAerror('CORA:wrongInputInConstructor',...
                        'Input "num" has a wrong format!'));
            end
        end

        % construct geometric set
        if isa(num, 'atomicProposition')
            obj.lhs = num;
            obj.variables = {name};
            obj.var = name;
            obj.logic = true;

            return;
        end

        % construct variable names
        if ~isempty(num)
        
            obj.variables = cell(num,1);

            for i = 1:num
               obj.variables{i} = [name, num2str(i)];
            end
        else
            obj.variables = {name}; 
        end
    end

    % operator overloading
    function res = mpower(obj1,obj2)
        res = power(obj1,obj2);
    end

    function res = mtimes(obj1,obj2)
        res = times(obj1,obj2);
    end

    function res = rdivide(obj1,obj2)
        res = obj1 * obj2^(-1);
    end

    function res = mrdivide(obj1,obj2)
        res = obj1 * obj2^(-1);
    end

    % mathematical functions
    function res = sin(obj)
        res = nonlinearOperator(obj,'sin');
    end

    function res = cos(obj)
        res = nonlinearOperator(obj,'cos');
    end

    function res = tan(obj)
        res = nonlinearOperator(obj,'tan');
    end

    function res = exp(obj)
        res = nonlinearOperator(obj,'exp');
    end
    
    function res = sqrt(obj)
        res = nonlinearOperator(obj,'sqrt');
    end

    function res = log(obj)
        res = nonlinearOperator(obj,'log');
    end

    function res = asin(obj)
        res = nonlinearOperator(obj,'asin');
    end

    function res = acos(obj)
        res = nonlinearOperator(obj,'acos');
    end

    function res = atan(obj)
        res = nonlinearOperator(obj,'atan');
    end

    function res = sinh(obj)
        res = nonlinearOperator(obj,'sinh');
    end

    function res = cosh(obj)
        res = nonlinearOperator(obj,'cosh');
    end

    function res = tanh(obj)
        res = nonlinearOperator(obj,'tanh');
    end
end

methods (Static = true)
    Z = generateRandom(varargin)
    stl = bool2stl(bool)
end

end

% ------------------------------ END OF CODE ------------------------------
