classdef (InferiorClasses = {?interval}) affine < taylm
% affine arithmetic class.
%
% Syntax:
%       object constructor: obj = affine(int)
%       object constructor: obj = affine(int, name, opt_method, eps, tolerance)
%       object constructor: obj = affine(lower_b, upper_b)
%       object constructor: obj = affine(lower_b, upper_b, name, opt_method, eps, tolerance)
%
% Inputs:
%    int - an interval
%    name - a cell containing a name of a variable
%    lower_b - a double, the lower bound of an interval
%    upper_b - a double, the upper bound of an interval
%    opt_method - method used to calculate interval over-approximations of
%                 taylor models 
%                  'int': standard interval arithmetic (default)
%                  'bnb': branch and bound method is used to find min/max
%                  'bnbAdv': branch and bound with re-expansion of taylor models
%    eps - precision for the selected optimization method (opt_method = 'bnb', 
%          opt_method = 'bnbAdv' and opt_method = 'linQuad')
%    tolerance - monomials with coefficients smaller than this value are
%                moved to the remainder
%
% Outputs:
%    obj - Generated Object
%
% Examples: 
%    % create affine object and taylor model object
%    a = affine(interval(0,pi/2),'a','int',[],1e-8);
%    t = taylm(interval(0,pi/2),6,'a','int',[],1e-8);
%
%    % compare the results
%    int_a = interval(sin(a))
%    int_t = interval(sin(t))
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval

% Author:       Dmitry Grebenyuk, Niklas Kochdumper
% Written:      22-September-2017
% Last update:  08-April-2018 (NK, extended constructor syntax)
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)

end
    
methods
    %class constructor
    function obj = affine(varargin)
        
        % check user input
        if nargin < 1
            error('Wrong syntax. Type "help affine" for more information.')
        end
        
        % first input is an interval object
        if isa(varargin{1},'interval')
           
            int = varargin{1};
            parseIndex = 2;
            
        % interval is specified by its boundaries
        else
            
            if nargin < 2
                error('Wrong syntax. Type "help affine" for more information.')
            end
            
            int = interval(varargin{1},varargin{2});
            parseIndex = 3;
        end
        
        % default values for the object properties
        opt_method = 'int';
        eps = 0.001;
        tolerance = 1e-8;
        
        % generate variable names if they are not provided
        try
           if nargin < parseIndex
               names = genDefaultVarNames(int,[],inputname(1));
           else
               names = genDefaultVarNames(int,varargin{parseIndex},inputname(1));
           end
        catch ex
           error(ex.message);
        end
        
        % parse input arguments
        if nargin >= parseIndex + 1 && ~isempty(varargin{parseIndex+1})
            opt_method = varargin{parseIndex+1};
        end
        if nargin >= parseIndex + 2 && ~isempty(varargin{parseIndex+2})
            eps = varargin{parseIndex+2};
        end
        if nargin >= parseIndex + 3 && ~isempty(varargin{parseIndex+3})
            tolerance = varargin{parseIndex+3};
        end
        
        % check if the selected optimization mehtod is feasible
        % (optimization with 'linQuad' is not possible for class affine)
        if ~ischar(opt_method) || ~ismember(opt_method,{'int','bnb','bnbAdv'})
          error('Wrong value for input argument "opt_method"!'); 
        end
        
        % create the object by calling the constructor of the superclass
        % (taylor model with max_order = 1)      
        obj = obj@taylm(int, 1, names, opt_method, eps, tolerance);
    end
             
end
end

%------------- END OF CODE -------