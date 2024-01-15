classdef conHyperplane < contSet
% conHyperplane - object constructor for constrained hyperplanes
%
% Description:
%    This class represents constrained hyperplane objects defined as
%    {x | a*x = b, C*x <= d}.
%
% Syntax:
%    obj = conHyperplane(a,b)
%    obj = conHyperplane(a,b,C,d)
%
% Inputs:
%    a - normal vector of the hyperplane a*x = b
%    b - offset of the hyperplane a*x = b
%    C - constraint matrix for the inequality constraints C*x <= d
%    d - constraint vector for the inequality constraints C*x <= d
%
% Outputs:
%    obj - conHyperplane object
%
% Example:
%    a = [-0.5; -1; 0.1];
%    b = -1;
%    C = [-0.6 0.8 -1.7;...
%          0.6 0.5 -0.8];
%    d = [1; 0.5];
%    hyp = conHyperplane(a,b,C,d);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: halfspace, example_conHyperplane.m

% Authors:       Matthias Althoff, Niklas Kochdumper, Victor Gassmann
% Written:       10-August-2011
% Last update:   22-November-2019 (NK, renamed + added additional constructors)
%                02-May-2020 (MW, added property validation)
%                19-March-2021 (MW, error messages)
%                22-March-2021 (VG, added 1D case)
%                14-December-2022 (TL, property check in inputArgsCheck)
%                16-August-2023 (TL, removed 1D check)
%                23-December-2023 (MW, disable using halfspace object)
% Last revision: 16-June-2023 (MW, restructure using auxiliary functions)

% ------------------------------ BEGIN CODE -------------------------------


properties (SetAccess = private, GetAccess = public)
    a;  % equality constraint vector
    b;  % equality constraint offset
    C;  % inequality constraint matrix
    d;  % inequality constraint offset

    % legacy
    h;  % halfspace
end
    
methods
    
    % class constructor
    function obj = conHyperplane(varargin)

        % 0. avoid empty instantiation
        if nargin == 0
            throw(CORAerror('CORA:noInputInSetConstructor'));
        end
        
        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'conHyperplane')
            obj = varargin{1}; return
        end

        % 2. parse input arguments: varargin -> vars
        [a,b,C,d] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(a,b,C,d,nargin);

        % 4. compute properties
        [a,b,C,d] = aux_computeProperties(a,b,C,d);

        % 5. assign properties
        obj.a = a;
        obj.b = b;
        obj.C = C;
        obj.d = d;
        
    end
         
    % methods in seperate files
    n = dim(hyp)
    val = distance(hyp,S)
    res = isequal(hyp1,hyp2,varargin)
    res = isHyperplane(hyp)
    P = polytope(hyp)
    han = plot(hyp,varargin)
    hyp = lift_(hyp,N,dims)
    Sproj = projectOnHyperplane(hyp,S)
    [res,S] = representsa_(hyp,type,tol,varargin)
        
    % display functions
    display(hyp)

end

methods (Static = true)
    hyp = generateRandom(varargin) % generate random conHyperplane
    hyp = empty(n) % instantiates an empty constrained hyperplane
    hyp = Inf(n) % instantiates a fullspace constrained hyperplane
end

% getter & setter ---------------------------------------------------------

methods
    % getter & setter for legacy Grest property
    function h = get.h(obj)
        warning(['CORA: The property conHyperplane.h is deprecated (since CORA 2024.1.0) and will be removed in a future release. ' ...
            'Please use conHyperplane.a and conHyperplane.b instead. ' ...
            'This change was made to be avoid confusion between the inequality constraint defining a halfspace and the equality constraint defining a constrained hyperplane.']);

        h = halfspace(obj.a,obj.b);
    end

    function obj = set.h(obj, h)
        warning(['CORA: The property polyZonotope.Grest is deprecated (since CORA 2024.1.0) and will be removed in a future release. ' ...
            'Please use conHyperplane.a and conHyperplane.b instead. ' ...
            'This change was made to be avoid confusion between the inequality constraint defining a halfspace and the equality constraint defining a constrained hyperplane.']);

        obj.a = h.c;
        obj.b = h.d;
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [a,b,C,d] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % check number of input arguments
    if nargin > 4
        throw(CORAerror('CORA:tooManyInputArgs',4));
    end

    % set default values depending on nargin
    [a,b,C,d] = setDefaultValues({[],[],[],[]},varargin);

end

function aux_checkInputArgs(a,b,C,d,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0
            
        inputArgsCheck({ ...
            {a, 'att', 'numeric', 'finite'}; ...
            {b, 'att', 'numeric', 'finite'}; ...
            {C, 'att', 'numeric', {'finite', 'matrix'}}; ...
            {d, 'att', 'numeric', 'finite'}; ...
        })

        if ~isvector(a) || ~isscalar(b)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Only a single equality constraint ax = b is supported.'));
        end

        if ~isempty(C)
            if length(a) ~= size(C,2)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    ['The dimension of the constraint matrix does not '...
                    'match the dimension of the halfspace.']));
            elseif size(C,1) ~= length(d)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    ['The length of the constraint offset does not '...
                    'match the dimension of the constraint matrix.']));
            end
        end
    end
    
end

function [a,b,C,d] = aux_computeProperties(a,b,C,d)

    % a must be a row vector
    a = reshape(a,1,[]);

    % d must be a column vector
    d = reshape(d,[],1);

end

% ------------------------------ END OF CODE ------------------------------
