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
%               14-December-2022 (TL, property check in inputArgsCheck)
% Last revision:---

%------------- BEGIN CODE --------------


properties (SetAccess = private, GetAccess = public)
    h; % halfspace obj
    C;
    d;
end
    
methods
    
    % class constructor
    function obj = conHyperplane(varargin)
        
        % parse input
        if nargin > 4
            % too many input arguments
            throw(CORAerror('CORA:tooManyInputArgs',4));
        end

        % default values
        h = []; C = []; d = 0;
        
        if nargin == 0
            h = halfspace();
            
        elseif nargin == 1 || nargin == 3
            if isa(varargin{1},'conHyperplane')
                % copy constructor
                obj = varargin{1};
                return
            end
            h = varargin{1};
            varargin = varargin(2:end);
            
        elseif nargin == 2 || nargin == 4
            h = halfspace(varargin{1},varargin{2});
            varargin = varargin(3:end);
        end

        % set default values
        [C,d] = setDefaultValues({C,d}, varargin);
        
        % check input arguments
        inputArgsCheck({ ...
            {h, 'att', 'halfspace'}; ...
            {C, 'att', 'numeric', {'finite', 'matrix'}}; ...
            {d, 'att', 'numeric', {'finite', 'column'}}; ...
        })

        if ~isempty(C)
            if length(h.c) ~= size(C,2)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'The dimension of the constraint matrix does not match the dimension of the halfspace.'));
            elseif size(C,1) ~= length(d)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'The length of the constraint offset does not match the dimension of the constraint matrix.'));
            end
        end
        
        % handle 1D case
        if dim(h)==1 && h.c~=0 && ~isempty(C)
            x = h.d/h.c;
            % check if x is in {x|Cx\leq d}
            X = mptPolytope(C,d);
            if ~contains_(X,x,'exact',eps)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    ['Assignment not consistent: implicit value for x ' ...
                    'given by hyperplane not contained in {x | C*x <= d}!']));
            end         
        end

        % set properties
        obj.h = h;
        obj.C = C;
        obj.d = d;
        
    end
         
    % methods in seperate files
    n = dim(hyp)
    val = distance(hyp,S)
    res = isempty(hyp)
    res = isequal(hyp1,hyp2,varargin)
    res = isHyperplane(hyp)
    P = mptPolytope(hyp)
    han = plot(hyp,varargin)
    hyp = projectHighDim(hyp,N,dims)
    Sproj = projectOnHyperplane(hyp,S)
        
    % display functions
    display(hyp)

end
end

%------------- END OF CODE --------------
