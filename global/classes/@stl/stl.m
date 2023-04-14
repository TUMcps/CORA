classdef stl
% stl - class representing a Signal Temporal Logic (STL) formula
%
% Syntax:
%    obj = stl(name)
%    obj = stl(name,num)
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

% Author:       Niklas Kochdumper
% Written:      09-November-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

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
end
    
methods
    
    % class constructor
    function obj = stl(name,varargin)
        
        % parse and check input arguments
        if ~ischar(name)
           throw(CORAerror('CORA:wrongInputInConstructor',...
                        'Input "name" has to be a string!'));
        end
        
        num = [];
        
        if nargin > 1
            
            num = varargin{1};
        
            if ~isnumeric(num) || ~isscalar(num) || floor(num) ~= ceil(num)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                        'Input "num" has a wrong format!'));
            end
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
end

methods (Static = true)
    Z = generateRandom(varargin)
    stl = bool2stl(bool)
end

end

%------------- END OF CODE --------------