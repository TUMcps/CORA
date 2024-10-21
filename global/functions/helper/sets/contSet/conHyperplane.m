function P = conHyperplane(varargin)
% conHyperplane - (deprecated) object constructor for constrained hyperplanes
%
% Description:
%    This class represents constrained hyperplane objects defined as
%    {x | a*x = b, C*x <= d}.
%
% Syntax:
%    P = conHyperplane(c,d,A,b)
%
% Inputs:
%    a - normal vector of the hyperplane a*x = b
%    b - offset of the hyperplane a*x = b
%    C - constraint matrix for the inequality constraints C*x <= d
%    d - constraint vector for the inequality constraints C*x <= d
%
% Outputs:
%    P - polytope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polytope

% Authors:       Tobias Ladner
% Written:       09-October-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

CORAwarning('CORA:deprecated','class','conHyperplane','CORA v2025', ...
    'When updating the code, please replace every function call ''conHyperplane(c,d,A,b)'' with ''polytope(A,b,c,d)''.', ...
    ['Constrained hyperplanes are a special case of polytopes. ' ...
    'As the benefit of having an additional class for this special case is minor,\n' ...
    'we removed it to improve maintainability. ' ...
    'The object is returned as a polytope.'] ...
);

% set default values
[a,b,C,d] = setDefaultValues({[],[],[],[]},varargin);

% make sure 'a' is a row vector
a = reshape(a,1,[]);

% call polytope constructor
P = polytope(C,d,a,b);

end

% ------------------------------ END OF CODE ------------------------------
