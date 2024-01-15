classdef zoo
% zoo - range bounding using intervals and Taylor models in parallel for
%    higher precision
%
% Syntax:
%    obj = zoo(int, methods)
%    obj = zoo(int, methods, names, max_order, eps, tolerance)
%
% Inputs:
%    int - interval object 
%    methods - cell-array containing the methods used in parallel as
%              strings (possible values: 'taylm(int)', 'taylm(bnb)',
%                       'taylm(bnbAdv)', 'taylm(linQuad)', 
%                       'affine(int)', 'affine(bnb)', 'affine(bnbAdv)', 
%                       and 'interval')
%    max_order - the maximal order of a polynomial stored in a polynomial part
%    eps - precision for the branch and bound (opt_method = 'bnb')
%          optimization
%    tolerance - monomials with coefficients smaller than this value are
%                moved to the remainder
%
% Outputs:
%    obj - generated object
%
% Examples:
%    % compute function bounds with interval arithmetic
%    I = interval([0;1],[2;4]);
%    int_i = I(1) * (I(1)-I(2)) + I(1)*I(2)
% 
%    % compute function bounds with class zoo
%    methods = {'taylm(int)'; 'affine(bnb)'; 'interval'};
%    z = zoo(I,methods,{'x';'y'},6,0.001,1e-8);
%    temp = z(1) * (z(1)-z(2)) + z(1)*z(2);
%    int_zoo = interval(temp)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, taylm

% Authors:       Dmitry Grebenyuk, Niklas Kochdumper
% Written:       05-November-2017
% Last update:   10-April-2018 (NK, modified object properties)
% Last revision: 16-June-2023 (MW, restructure using auxiliary functions)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    
    % cell array containing the names of the applied methods as strings
    method;
    
    % cell array containing the class objects for the applied methods
    objects;
    
end
    
methods
    
    % class constructor
    function obj = zoo(varargin)

        % 0. avoid empty instantiation
        if nargin == 0
            throw(CORAerror('CORA:noInputInSetConstructor'));
        end

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'zoo')
            obj = varargin{1}; return
        end

        % 2. parse input arguments: varargin -> vars
        [int,methods,names,max_order,eps,tolerance] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(int,methods,names,max_order,eps,tolerance,nargin);

        % 4. compute object
        [method,objects] = aux_computeObject(int,methods,names,max_order,eps,tolerance,inputname(1));

        % 5. assign properties
        obj.method = method;
        obj.objects = objects;

    end
    
    % class methods
    function res = acos(obj); res = zooComputation(@acos,obj); end
    function res = asin(obj); res = zooComputation(@asin,obj); end
    function res = atan(obj); res = zooComputation(@atan,obj); end
    function res = cos(obj); res = zooComputation(@cos,obj); end
    function res = cosh(obj); res = zooComputation(@cosh,obj); end
    function res = exp(obj); res = zooComputation(@exp,obj); end
    function res = log(obj); res = zooComputation(@log,obj); end
    function res = sin(obj); res = zooComputation(@sin,obj); end
    function res = sinh(obj); res = zooComputation(@sinh,obj); end
    function res = sqrt(obj); res = zooComputation(@sqrt,obj); end
    function res = tan(obj); res = zooComputation(@tan,obj); end
    function res = tanh(obj); res = zooComputation(@tanh,obj); end
    function res = mpower(obj1,obj2); res = zooComputation(@mpower,obj1,obj2); end
    function res = power(obj1,obj2); res = zooComputation(@power,obj1,obj2); end

    function res = isemptyobject(obj); res = false; end
             
end

end


% Auxiliary functions -----------------------------------------------------

function [int,method,names,max_order,eps,tolerance] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % check number of input arguments
    if nargin > 6
        throw(CORAerror('CORA:tooManyInputArgs',6));
    end
    if nargin < 2
        throw(CORAerror('CORA:notEnoughInputArgs',2));
    end

    % set default values (first two always given)
    [int,method,names,max_order,eps,tolerance] = ...
        setDefaultValues({[],[],[],6,0.001,1e-8},varargin);

end

function aux_checkInputArgs(int,methods,names,max_order,eps,tolerance,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0

        % int has to be an interval
        if ~isa(int,'interval')
            throw(CORAerror('CORA:wrongValue','first',...
                'has to be an interval object.'));
        end

        % check methods
        if ~all(cell2mat(cellfun(@ischar,methods,'UniformOutput',false)))
            throw(CORAerror('CORA:wrongValue','second',"has to be a cell-array of methods."));
        end
        
        temp = ismember(methods,{'taylm(int)','taylm(bnb)','taylm(bnbAdv)', ...
                'taylm(linQuad)','affine(int)','affine(bnb)','affine(bnbAdv)', ...
                'interval'});
        if ~all(temp)
           ind = find(temp == 0);
           str = methods{ind(1)};
           throw(CORAerror('CORA:wrongValue','second',...
               "'taylm(int)', 'taylm(bnb)', 'taylm(bnbAdv)', 'taylm(linQuad)', " + ...
               "'affine(int)', 'affine(bnb)', 'affine(bnbAdv)', 'interval'"));
        end

        % correct value for max_order
        if ~isnumeric(max_order) || ~isscalar(max_order) || mod(max_order,1) ~= 0
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Maximum order must be an integer greater than zero.'));
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

function [method,objects] = aux_computeObject(int,methods,names,max_order,eps,tolerance,varname)

    % generate variable names if they are not provided
    names = genDefaultVarNames(int,names,varname);

    % sort methods alphabetically
    method = sort(methods);
    % generate the objects
    objects = cell(length(methods),1);
    
    for i = 1:length(method)
       
        m = method{i};
        
        switch m
            
            case 'taylm(int)'
                objects{i} = taylm(int,max_order,names,'int',eps,tolerance);
            
            case 'taylm(bnb)'
                objects{i} = taylm(int,max_order,names,'bnb',eps,tolerance);
                
            case 'taylm(bnbAdv)'
                objects{i} = taylm(int,max_order,names,'bnbAdv',eps,tolerance);
            
            case 'taylm(linQuad)'
                objects{i} = taylm(int,max_order,names,'linQuad',eps,tolerance);
                
            case 'affine(int)'
                objects{i} = affine(int,names,'int',eps,tolerance);
                
            case 'affine(bnb)'
                objects{i} = affine(int,names,'bnb',eps,tolerance);
                
            case 'affine(bnbAdv)'
                objects{i} = affine(int,names,'bnbAdv',eps,tolerance);
                
            case 'interval'
                objects{i} = int;
            
        end
    end

end

% ------------------------------ END OF CODE ------------------------------
