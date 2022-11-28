classdef conHyperplane < contSet
% conHyperplane - object constructor for constrained hyperplanes
%
% Description:
%    This class represents constrained hyperplane objects defined as
%    {x | a*x = b, C*x <= d}.
%
% Syntax:  
%    obj = conHyperplane(hs)
%    obj = conHyperplane(a,b)
%    obj = conHyperplane(hs,C,d)
%    obj = conHyperplane(a,b,C,d)
%
% Inputs:
%    hs - halfspace object defining the constraint a*x = b
%    a - normal vector of the hyperplane a*x = b
%    b - offset of the hyperplane a*x = b
%    C - constraint matrix for the inequality constraints C*x <= d
%    d - constraint vector for the inequality constraints C*x <= d
%
% Outputs:
%    obj - conHyperplane object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: halfspace, example_conHyperplane.m

% Author:       Matthias Althoff, Niklas Kochdumper, Victor Gassmann
% Written:      10-August-2011
% Last update:  22-Nov-2019 (NK, renamed + added additional constructors)
%               02-May-2020 (MW, added property validation)
%               19-March-2021 (MW, error messages)
%               22-March-2021 (VG, added 1D case)
% Last revision:---

%------------- BEGIN CODE --------------


properties (SetAccess = private, GetAccess = public)
    h = []; % halfspace obj
    C (:,:) {mustBeNumeric,mustBeFinite} = [];
    d (:,1) {mustBeNumeric,mustBeFinite} = 0;
end
    
methods
    
    % class constructor
    function obj = conHyperplane(varargin)
        
        if nargin == 0
            obj.h = halfspace();
        
        elseif nargin == 1
            % copy constructor
            if isa(varargin{1},'conHyperplane')
                obj = varargin{1};
            elseif isa(varargin{1},'halfspace')
                obj.h = varargin{1};
            else
                throw(CORAerror('CORA:wrongValue','first',...
                    "'conHyperplane' object or 'halfspace' object"));
            end
            
        elseif nargin == 2
            obj.h = halfspace(varargin{1},varargin{2});
            
        elseif nargin == 3
            if ~isa(varargin{1},'halfspace') 
                throw(CORAerror('CORA:wrongInputInConstructor',...
                     'If three input arguments are given, the first one has to be a ''halfspace'' object'));
            elseif dim(varargin{1}) ~= size(varargin{2},2)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'The dimension of the constraint matrix does not match the dimension of the halfspace.'));
            elseif size(varargin{2},1) ~= length(varargin{3})
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'The length of the constraint offset does not match the dimension of the constraint matrix.'));
            else
                obj.h = varargin{1};
                obj.C = varargin{2};
                obj.d = varargin{3};
            end
            
        elseif nargin == 4
            if length(varargin{1}) ~= size(varargin{3},2)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'The dimension of the constraint matrix does not match the dimension of the halfspace.'));
            elseif size(varargin{3},1) ~= length(varargin{4})
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'The length of the constraint offset does not match the dimension of the constraint matrix.'));
            else
                obj.h = halfspace(varargin{1},varargin{2});
                obj.C = varargin{3};
                obj.d = varargin{4};
            end
            
        elseif nargin > 4
            % too many input arguments
            throw(CORAerror('CORA:tooManyInputArgs',4));
        end
        
        % handle 1D case
        if dim(obj)==1 && obj.h.c~=0 && ~isempty(obj.C)
            x = obj.h.d/obj.h.c;
            % check if x is in {x|Cx\leq d}
            X = mptPolytope(obj.C,obj.d);
            if ~contains(X,x)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Assignment not consistent: implicit value for x given by hyperplane not contained in {x | C*x <= d}!'));
            end         
        end
        
        % set parent object properties
        obj.dimension = dim(obj.h);
        
    end
         
    % methods in seperate files  
    res = and(hyp,S)
    res = contains(hyp,S)
    n = dim(hyp)
    val = distance(hyp,S)
    res = isempty(hyp)
    res = isequal(hyp1,hyp2,varargin)
    res = isHyperplane(hyp)
    res = isIntersecting(hyp,S,varargin)
    P = mptPolytope(hyp)
    han = plot(hyp,varargin)
    hyp = projectHighDim(hyp,N,dims)
    Sproj = projectOnHyperplane(hyp,S)
    [val,x] = supportFunc(hyp,d,type)
        
    % display functions
    display(hyp)

end
end

%------------- END OF CODE --------------
