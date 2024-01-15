classdef (InferiorClasses = {?interval}) affine < taylm
% affine arithmetic class.
%
% Syntax:
%    obj = affine(lb, ub)
%    obj = affine(lb, ub, name, opt_method, eps, tolerance)
%
% Inputs:
%    name - a cell containing a name of a variable
%    lb - lower bound of an interval
%    ub - upper bound of an interval
%    opt_method - method used to calculate interval over-approximations of
%                 taylor models 
%                  'int': standard interval arithmetic (default)
%                  'bnb': branch and bound method is used to find min/max
%                  'bnbAdv': branch and bound with re-expansion of Taylor models
%    eps - precision for the selected optimization method (opt_method = 'bnb', 
%          opt_method = 'bnbAdv', and opt_method = 'linQuad')
%    tolerance - monomials with coefficients smaller than this value are
%                moved to the remainder
%
% Outputs:
%    obj - Generated Object
%
% Examples:
%    % create affine object and taylor model object
%    a = affine(0,pi/2,'a','int',[],1e-8);
%    t = taylm(interval(0,pi/2),6,'a','int',[],1e-8);
%
%    % compare the results
%    int_a = interval(sin(a))
%    int_t = interval(sin(t))
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval

% Authors:       Dmitry Grebenyuk, Niklas Kochdumper
% Written:       22-September-2017
% Last update:   08-April-2018 (NK, extended constructor syntax)
%                24-July-2023 (MW, integrate isemptyobject_)
% Last revision: 16-June-2023 (MW, restructure using auxiliary functions)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    % inherits properties from taylm
end
    
methods
    %class constructor
    function obj = affine(varargin)

        % 0. avoid empty instantiation
        if nargin == 0
            throw(CORAerror('CORA:noInputInSetConstructor'));
        end

        % 1. copy constructor: not allowed due to obj@taylm below
%         if nargin == 1 && isa(varargin{1},'affine')
%             obj = varargin{1}; return
%         end

        % 1. parse input arguments: varargin -> vars
        [int,name,opt_method,eps,tolerance] = aux_parseInputArgs(varargin{:});

        % 2. check correctness of input arguments
        aux_checkInputArgs(int,name,opt_method,eps,tolerance,nargin);

        % 3. compute properties
        [int,name,opt_method,eps,tolerance] = ...
            aux_computeProperties(int,name,opt_method,eps,tolerance,inputname(1));
        
        % 4. assign properties
        % ...create the object by calling the constructor of the superclass
        % (taylor model with max_order = 1)
        obj = obj@taylm(int,1,name,opt_method,eps,tolerance);
    end

    
    % conversion of affine object to taylm object
    function tay = taylm(aff)
        
        % copy all properties
        c = metaclass(aff);
        prop = c.Properties;
        
        for i = 1:length(prop)
            if ~prop{i}.Dependent
                tay.(prop{i}.Name) = aff.(prop{i}.Name);
            end
        end

    end

    % determine whether an affine object is fully empty
    function res = isemptyobject(aff)
        % no empty object allowed
        res = false;
    end
             
end
end


% Auxiliary functions -----------------------------------------------------

function [int,name,opt_method,eps,tolerance] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % check number of input arguments
    if nargin > 6
        throw(CORAerror('CORA:tooManyInputArgs',6));
    end

    % no input arguments
    if nargin == 0
        int = []; name = []; opt_method = []; eps = []; tolerance = [];
        return
    end

    % set default values
%     if isa(varargin{1},'interval')
%         [int,name,opt_method,eps,tolerance] = ...
%             setDefaultValues({[],[],'int',0.001,1e-8},varargin);
%     else
        [lb,ub,name,opt_method,eps,tolerance] = ...
            setDefaultValues({[],[],[],'int',0.001,1e-8},varargin);
        % init interval from lower and upper bounds
        int = interval(lb,ub);
%     end

end

function aux_checkInputArgs(int,name,opt_method,eps,tolerance,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0

        % check if the selected optimization method is feasible
        % (optimization with 'linQuad' is not possible for class affine)
        if ~ischar(opt_method) || ~ismember(opt_method,{'int','bnb','bnbAdv'})
            throw(CORAerror('CORA:wrongValue','third/fourth',"'int','bnb' or 'bnbAdv'"));
        end
        
    end

end

function [int,name,opt_method,eps,tolerance] = ...
    aux_computeProperties(int,name,opt_method,eps,tolerance,varname)
% compute properties of affine object

    % generate variable names
    name = genDefaultVarNames(int,name,varname);

end

% ------------------------------ END OF CODE ------------------------------
