classdef conHyperplane
% conHyperplane class 
%
% Description:
%    This class respresnts constrained hyperplane objects defined as
%    {x | a*x = b, C*x <= d}
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
%    C - constrained matrix for the inequality constraints C*x <= d
%    d - constrained vector for the inequality constraints C*x <= d
%
% Outputs:
%    obj - generated conHyperplane object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: halfspace

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      10-August-2011
% Last update:  22-Nov-2019 (NK, renamed + added additional constructors)
%               02-May-2020 (MW, added property validation)
% Last revision:---

%------------- BEGIN CODE --------------


properties (SetAccess = private, GetAccess = public)
    h = []; %halfspace obj
    C (:,:) {mustBeNumeric,mustBeFinite} = [];
    d (:,1) {mustBeNumeric,mustBeFinite} = 0;
end
    
methods
    
    % class constructor
    function obj = conHyperplane(varargin)
        
        if nargin == 1
            % copy constructor
            if isa(varargin{1},'conHyperplane')
                obj = varargin{1};
            else
                obj.h = varargin{1};
            end
            
        elseif nargin == 2
            obj.h = halfspace(varargin{1},varargin{2});
            
        elseif nargin == 3
            obj.h = varargin{1};
            obj.C = varargin{2};
            obj.d = varargin{3};
            
        elseif nargin == 4
            obj.h = halfspace(varargin{1},varargin{2});
            obj.C = varargin{3};
            obj.d = varargin{4};
        end
        
    end
         
    % methods in seperate files  
    obj = and(obj,S)
    res = isempty(hyp)
    res = isequal(hyp1,hyp2)
    res = isIntersecting(obj1,obj2,varargin)
    P = mptPolytope(obj)
    han = plot(obj,varargin)
    res = projectHighDim(obj,N,dims)
    res = projectOnHyperplane(h, S)    
        
    % display functions
    display(obj)

end
end

%------------- END OF CODE --------------
