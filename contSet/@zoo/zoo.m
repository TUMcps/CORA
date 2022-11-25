classdef zoo
% Zoo class computes intervals and Taylor models in parallel for better
% precision
%
% Syntax:
%       object constructor: obj = zoo(int, methods)
%                           obj = zoo(int, methods, names, max_order, eps, tolerance)
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
%
% Outputs:
%    obj - Generated Object
%
% Examples:
%   % compute function bounds with interval arithmetic
%   i = interval([0;1],[2;4]);
%   int_i = i(1) * (i(1)-i(2)) + i(1)*i(2)
%
%   % compute function bounds with class zoo
%   methods = {'taylm(int)'; 'affine(bnb)'; 'interval'};
%   z = zoo(i,methods,{'x';'y'},6,0.001,1e-8);
%   int_zoo = interval(z(1) * (z(1)-z(2)) + z(1)*z(2))
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, taylm

% Author:       Dmitry Grebenyuk, Niklas Kochdumper
% Written:      05-November-2017
%
% Last update:  10-April-2018 (NK, modified object properties)
%               ---
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    
    % cell array containing the names of the applied methods as strings
    method
    
    % cell array contaiing the class objects for the applied methods
    objects
    
end
    
methods
    
    % class constructor
    function obj = zoo(int, methods, varargin)
        
        % default settings
        max_order = 6;
        tolerance = 1e-8;
        eps = 0.001;
        
        % parse input arguments
        if nargin < 2 || nargin > 6
           error('Wrong syntax! Type "help zoo" for more information.'); 
        end
        
        if nargin >= 4 && ~isempty(varargin{2})
           max_order = varargin{2}; 
        end
        if nargin >= 5 && ~isempty(varargin{3})
           eps = varargin{3}; 
        end
        if nargin >= 6 && ~isempty(varargin{4})
           tolerance = varargin{4}; 
        end
        
        % generate variable names if they are not provided
        try
           if nargin < 3
               names = genDefaultVarNames(int,[],inputname(1));
           else
               names = genDefaultVarNames(int,varargin{1},inputname(1));
           end
        catch ex
           error(ex.message);
        end
        
        
        % check user input
        if ~all(cell2mat(cellfun(@ischar,methods,'UniformOutput',false))) || ~isa(int,'interval')
            error('Wrong syntax! Type "help zoo" for more information.');
        else
            temp = ismember(methods,{'taylm(int)','taylm(bnb)','taylm(bnbAdv)', ...
                    'taylm(linQuad)','affine(int)','affine(bnb)','affine(bnbAdv)', ...
                    'interval'});
            if ~all(temp)
               ind = find(temp == 0);
               str = methods{ind(1)};
               error('The string %s is not a valid value for input argument name!',str);
            end 
        end
        
        % generate the objects
        obj.method = sort(methods);     % sort alphabetically
        obj.objects = cell(length(methods),1);
        
        for i = 1:length(obj.method)
           
            m = obj.method{i};
            
            switch m
                
                case 'taylm(int)'
                    obj.objects{i} = taylm(int,max_order,names,'int',eps,tolerance);
                
                case 'taylm(bnb)'
                    obj.objects{i} = taylm(int,max_order,names,'bnb',eps,tolerance);
                    
                case 'taylm(bnbAdv)'
                    obj.objects{i} = taylm(int,max_order,names,'bnbAdv',eps,tolerance);
                
                case 'taylm(linQuad)'
                    obj.objects{i} = taylm(int,max_order,names,'linQuad',eps,tolerance);
                    
                case 'affine(int)'
                    obj.objects{i} = affine(int,names,'int',eps,tolerance);
                    
                case 'affine(bnb)'
                    obj.objects{i} = affine(int,names,'bnb',eps,tolerance);
                    
                case 'affine(bnbAdv)'
                    obj.objects{i} = affine(int,names,'bnbAdv',eps,tolerance);
                    
                case 'interval'
                    obj.objects{i} = int;
                
            end
        end

    end
    
    % class methods
    function res = acos(obj); res = zooComputation(@acos, obj); end
    function res = asin(obj); res = zooComputation(@asin, obj); end
    function res = atan(obj); res = zooComputation(@atan, obj); end
    function res = cos(obj); res = zooComputation(@cos, obj); end
    function res = cosh(obj); res = zooComputation(@cosh, obj); end
    function res = exp(obj); res = zooComputation(@exp, obj); end
    function res = log(obj); res = zooComputation(@log, obj); end
    function res = sin(obj); res = zooComputation(@sin, obj); end
    function res = sinh(obj); res = zooComputation(@sinh, obj); end
    function res = sqrt(obj); res = zooComputation(@sqrt, obj); end
    function res = tan(obj); res = zooComputation(@tan, obj); end
    function res = tanh(obj); res = zooComputation(@tanh, obj); end
    function res = mpower(obj1, obj2); res = zooComputation(@mpower, obj1, obj2); end
    function res = power(obj1, obj2); res = zooComputation(@power, obj1, obj2); end
             
end
end

%------------- END OF CODE -------