function P = halfspace(varargin)
% halfspace - (deprecated) object constructor for halfspaces 
%
% Description:
%    This class represents halfspace objects defined as
%    {x | c'*x <= d}.
%
% Syntax:
%    P = halfspace(c,d)
%
% Inputs:
%    c - normal vector
%    d - offset
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

CORAwarning('CORA:deprecated','class','halfspace','CORA v2025', ...
    'When updating the code, please replace every function call ''halfspace(c,d)'' with ''polytope(c,d)''.', ...
    ['Halfspaces are a special case of polytopes. ' ...
    'As the benefit of having an additional class for this special case is minor,\n' ...
    'we removed it to improve maintainability. ' ...
    'The object is returned as a polytope.'] ...
);

% set default values
[c,d] = setDefaultValues({[],[]},varargin);

% make sure 'a' is a row vector
c = reshape(c,1,[]);

% call polytope constructor
P = polytope(c,d);

end

% ------------------------------ END OF CODE ------------------------------
