classdef (InferiorClasses = {?interval}) taylm
% taylm (Taylor model) class.
%
% Syntax:  
%    obj = taylm()
%    obj = taylm(int)
%    obj = taylm(int,max_order,names,opt_method,eps,tolerance)
%    obj = taylm(func,int)
%    obj = taylm(func,int,max_order,opt_method,eps,tolerance)
%
% Inputs:
%    int - interval object that defines the ranges of the variables
%    max_order - the maximal order of a polynomial stored in a polynomial part
%    names - cell-array containing the names of the symbolic variables as
%            characters (same size as interval matrix 'int')
%    func - symbolic function 
%    opt_method - method used to calculate interval over-approximations of
%                 taylor models 
%                  'int': standard interval arithmetic (default)
%                  'bnb': branch and bound method is used to find min/max
%                  'bnbAdv': branch and bound with re-expansion of taylor models
%                  'linQuad': optimization with Linear Dominated Bounder (LDB)
%                             and Quadratic Fast Bounder (QFB)
%    eps - precision for the selected optimization method (opt_method = 'bnb', 
%          opt_method = 'bnbAdv' and opt_method = 'linQuad')
%    tolerance - monomials with coefficients smaller than this value are
%                moved to the remainder
%
% Outputs:
%    obj - generated object
%
% Examples: 
%    % create and manipulate simple taylor models
%    tx = taylm(interval(1,2),4,'x');
%    ty = taylm(interval(3,4),4,'y');
%    t = sin(tx+ty) + exp(-tx) + ty*tx;
%    interval(t)
%
%    % create a vector of taylor models
%    tvec = taylm(interval([1;3],[2;4]),4,{'x';'y'});
%    t = sin(tvec(1)+tvec(2)) + exp(-tvec(1)) + tvec(1)*tvec(2);
%    interval(t)
%
%    % create a taylor model from a symbolic function
%    syms x y
%    func = sin(x+y) + exp(-x) + x*y;
%    t = taylm(func,interval([1;3],[2;4]),4);
%    interval(t)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval
%
% References: 
%   [1] K. Makino et al. "Taylor Models and other validated functional 
%       inclusion methods"
%   [2] M. Althoff et al. "Implementation of Taylor models in CORA 2018
%       (Tool Presentation)"

% Author:       Dmitry Grebenyuk, Niklas Kochdumper
% Written:      29-March-2016             
% Last update:  18-July-2017 (DG) Multivariable polynomial pack is added
%               29-July-2017 (DG, NK) The NK' code is imerged with the DG'
%               11-October-2017 (DG) Syms as an input
%               3-April-2018 (NK) Restructured constructor
%               02-May-2020 (MW) add property validation 
% Last revision:16-June-2023 (MW, restructure using auxiliary functions)

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    % coefficients of polynomial terms. 
    % format:       column vector
    coefficients = [];
    
    % monomials for of poly. terms (i.e. x*(y.^2)*z -> [1;2;1])
    % format:      cell array of column vectors
    monomials = [];
    
    % remainder term of the Taylor model
    % format:       interval
    remainder = interval(0,0);
    
    % list with names of the symbolic variables, i.e. {'x','y','z'} or
    % {'x'}. When created, it has one name, the name of int.
    % format:       cell array of strings
    names_of_var = {};
    
    % stores the maximal order of a polynomial
    % format:       integer
    max_order = 6; 
    
    % defines the method used to determine an interval over-approximation
    % format:       string ('int', 'bnb', 'bnbAdv', 'linQuad')
    opt_method = 'int';
    
    % precision for the branch and bound optimization (opt_method = 'bnb')
    % format:       scalar (> 0)
    eps = 0.001;
    
    % coefficients smaller than this value get moved to the remainder
    % format:       scalar (> 0)
    tolerance = 1e-8;
    
end
    
methods
    % class constructor
    function obj = taylm(varargin)

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'taylm')
            obj = varargin{1}; return
        end

        % 2. parse input arguments: varargin -> vars
        [func,int,max_order,names,opt_method,eps,tolerance] = ...
            aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(func,int,max_order,names,opt_method,eps,tolerance,nargin);

        % 4. compute object
        obj = aux_computeObject(obj,func,int,max_order,names,opt_method,eps,tolerance,inputname(1));

    end
    

    % methods in seperate files 
    res = plus(summand1,summand2)
    res = minus(minuend,subtrahend)
    res = times(factor1, factor2)
    res = mtimes(factor1,factor2)
    res = mrdivide(numerator,denominator)
    res = mpower(base,exponent)
    res = power(base,exponent)
    obj = uplus(obj)   % overloaded '+' operator for single operand
    obj = uminus(obj)  % overloaded '-' operator for single operand    
    
    newObj = subsref(obj, S) % retrieves values from arrays
    obj = subsasgn(obj, S, value) % assigns values to arrays
%     res = length(obj) % returns the length of the array
%     varargout = size(obj) % returns size of object
    tm = setName( tm, name )   % changes the name of a variable in a Taylor model
    coef = getCoef( obj )   % returns coefficients
    rem = getRem( obj )     % returns remainder
    res = det( obj )        % determinant of Taylor model matrix
    res = log( obj )        % compute 'log' for Taylor models
    res = sqrt( obj )       % compute formula the square root
    res = exp( obj )        % exponential function
    res = sin( obj )        % sine function
    res = cos( obj )        % cosine function
    res = tan( obj )        % overloaded 'tan()' 
    res = sinh( obj )       % overloaded 'sinh()'
    res = cosh( obj )       % overloaded 'cosh()' 
    res = tanh( obj )       % overloaded 'tanh()'
    res = asin( obj )       % overloaded 'asin()'
    res = acos( obj )       % overloaded 'acos()'
    res = atan( obj )       % overloaded 'asin()'
    res = getSyms( obj )    % returns a polynomial in a sym form
    res = trace(obj)        % trace for TM matrices
    res = prod(obj,varargin)    % product of array elements
    res = sum(obj,varargin)     % sum of array elements
    
    %display functions
    display(obj)
end

end


% Auxiliary Functions -----------------------------------------------------

function [func,int,max_order,names,opt_method,eps,tolerance] = ...
    aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % due to different supported syntaxes, number of allowed input
    % arguments depends on input arguments (see below)

    % default values
    def_max_order = 6;
    def_opt_method = 'int';
    def_eps = 0.001;
    def_tolerance = 1e-8;
    % init values
    func = []; int = []; names = {};

    % no input arguments
    if nargin == 0
        max_order = def_max_order; opt_method = def_opt_method;
        eps = def_eps; tolerance = def_tolerance;
        return
    end

    % convert numeric to interval
    if isnumeric(varargin{1})
        varargin{1} = interval(varargin{1});
    end

    % two cases: first input argument either interval or symbolic object
    if isa(varargin{1},'interval')
        % syntaxes:
        % - obj = taylm(int)
        % - obj = taylm(int,max_order,names,opt_method,eps,tolerance)

        % check number of input arguments
        if nargin > 6
            throw(CORAerror('CORA:tooManyInputArgs',6));
        end
        
        % read out values
        [int,max_order,names,opt_method,eps,tolerance] = ...
            setDefaultValues({[],def_max_order,[],def_opt_method,def_eps,def_tolerance},varargin);
        
        % set func to empty for checking of input arguments and computation
        % of object properties
        func = [];
    
    elseif isa(varargin{1},'sym')
        % syntaxes:
        % - obj = taylm(func,int)
        % - obj = taylm(func,int,max_order,opt_method,eps,tolerance)

       % check number of input arguments
       if nargin < 2
           throw(CORAerror('CORA:notEnoughInputArgs',2));
       elseif nargin > 6
           throw(CORAerror('CORA:tooManyInputArgs',6));
       end

       % read out values
       [func,int,max_order,opt_method,eps,tolerance] = ...
           setDefaultValues({[],[],def_max_order,def_opt_method,def_eps,def_tolerance},varargin);

       % set names to empty (always computed automatically if func given)
       names = [];

    end

    % if values where there are default values are [], choose default value
    % note: find better solution for this
    if isempty(max_order);  max_order = def_max_order; end
    if isempty(opt_method); opt_method = def_opt_method; end
    if isempty(eps);        eps = def_eps; end
    if isempty(tolerance);  tolerance = def_tolerance; end

end

function aux_checkInputArgs(func,int,max_order,names,opt_method,eps,tolerance,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0

        % correct value for max_order
        if ~isnumeric(max_order) || ~isscalar(max_order) || mod(max_order,1) ~= 0
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Maximum order must be an integer greater than zero.'));
        end

        % correct method for optimization
        if ~(ischar(opt_method) || isstring(opt_method)) || ~ismember(char(opt_method),{'int','bnb','bnbAdv','linQuad'})
            throw(CORAerror('CORA:wrongValue','fourth',"'int', 'bnb', 'bnbAdv', or 'linQuad'"));
        end

        % correct value for eps
        if ~isnumeric(eps) || ~isscalar(eps) || eps <= 0
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Precision for branch and bound optimization must be a scalar greater than zero.'));
        end

        % correct value for tolerance
        if ~isnumeric(tolerance) || ~isscalar(tolerance) || tolerance <= 0
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Tolerance be a scalar greater than zero.'));
        end

    end

end

function obj = aux_computeObject(obj,func,int,max_order,names,opt_method,eps,tolerance,varname)
% compute properties of taylm object

    if isempty(func) && isempty(int)
        % immediate exit
        obj.coefficients = 0;
        obj.max_order = max_order;
        obj.opt_method = opt_method;
        obj.eps = eps;
        obj.tolerance = tolerance;
        return
    end

    if isempty(func) && ~isempty(int)
        % varargin{1} was an interval

        % generate variable names if they are not provided
        names = genDefaultVarNames(int,names,varname);
        
        % generate the taylor model
        if isscalar(int)
        
            obj.coefficients = [center(int); rad(int)];
            obj.monomials = hashFunction([0; 1]);
            obj.max_order = max_order;
            obj.opt_method = opt_method;
            obj.eps = eps;
            obj.tolerance = tolerance;
            obj.remainder = interval(0,0);
            if ~iscell(names)
                obj.names_of_var = {names};
            else
                obj.names_of_var = names;
            end
        
        else
        
            % use "repelem" instead of "arrayfunc" to initialize the
            % object-array, because only this way the initialization of
            % subclass arrays with the superclass constructor is
            % suported
            obj = repelem(obj,size(int,1),size(int,2));
            for i = 1:size(int,1)
                for j = 1:size(int,2)
                    obj(i,j).coefficients = [center(int(i,j)); rad(int(i,j))];
                    obj(i,j).monomials = hashFunction([0; 1]);
                    obj(i,j).max_order = max_order;
                    obj(i,j).opt_method = opt_method;
                    obj(i,j).eps = eps;
                    obj(i,j).tolerance = tolerance;
                    obj(i,j).remainder = interval(0,0);
                    obj(i,j).names_of_var = names(i,j);
                end
            end
        
        end


    elseif ~isempty(func)
        % varargin{1} was a symbolic function

        % scalar vs. matrix input
        if ~isscalar(func)
        
            % assign each symbolic variable the correct interval of
            % values
            v = symvar(func);
            obj = cell(size(func));
            
            for i = 1:size(obj,1)
                for j = 1:size(obj,2)
                
                    v_ = symvar(func(i,j));
                    int_ = interval(zeros(length(v_),1),zeros(length(v_),1));
                    counter = 1;
                    k = 1;
                    while k <= length(v) && counter <= length(v_)
                        if isequal(v(k),v_(counter))
                            int_(counter) = int(k);
                            counter = counter + 1;
                        end
                            k = k + 1;
                    end
                    
                    % call taylm constructor for a scalar function
                    obj{i,j} = taylm(func(i,j),int_,max_order,opt_method,eps,tolerance);
                
                end
            end
            
            % bring the object to the correct format
            A = cat(1, obj{:});
            obj = reshape(A, size(obj));          
            
        else
            
            % extract variable names
            v = symvar(func);
            names = cellfun(@(x) char(x),num2cell(transpose(v)),'UniformOutput',false);
            
            if length(names) > length(int)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                ['The length of the input argument ''int'' '...
                'has to be identical to the number of variables in the symbolic function!'])); 
            end
            
            % create taylor models for all variables
            for i = 1:length(v)
                str = sprintf('t%s = taylm(interval(%i,%i),%i,''%s'',''%s'',%e,%e);', ...
                names{i},infimum(int(i)), ...
                supremum(int(i)),max_order,names{i},opt_method,eps,tolerance);
                tay(i) = sym(sprintf('t%s',names{i}));
                
                eval(str); 
            end
            
            % evaluate the symbolic formula with taylor models
            if ~isempty(v)
                temp = subs(func,v,tay);
                str = sprintf('obj = %s;',char(temp));
                eval(str);
            else
                temp = eval(func);
                obj = taylm(interval(temp,temp),max_order,'const',opt_method,eps,tolerance);
            end
    
        end

    end

end

%------------- END OF CODE -------