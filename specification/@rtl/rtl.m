classdef rtl
% rtl - class representing a Reachset Temporal Logic (RTL) formula
%
% Syntax:
%    obj = rtl(eq)
%
% Inputs:
%    eq - STL formula (class stl)
%
% Outputs:
%    obj - generated rtl variable
%
% Example:
%    x = stl('x',2);
%    eq = rtl(x(1) <= 5) | rtl(x(2) <= 3)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: rtl

% Authors:       Niklas Kochdumper
% Written:       09-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    type = 'all'            % type of the operator ('next', '&', etc.)
    lhs = []                % object forming left-hand side of operator
    rhs = []                % object forming right-hand side of operator
    time = []               % time shift for next-opearator
end
    
methods
    
    % class constructor
    function obj = rtl(eq)
        
        % parse and check input arguments
        if ~isa(eq,'stl')
            throw(CORAerror('CORA:wrongInputInConstructor',...
                      'Input "eq" has to be an object of class "stl"!'));
        end
        
        obj.lhs = eq;
    end
end
end

% ------------------------------ END OF CODE ------------------------------
