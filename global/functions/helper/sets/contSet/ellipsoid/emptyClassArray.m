function A = emptyClassArray(classname,varargin)
% emptyClassArray - constructs class array of specified dimensions
%
% Syntax:
%    A = emptyClassArray(classname,n)
%    A = emptyClassArray(classname,n1,n2,...)
%    A = emptyClassArray(classname,[n1,n2,...])
%
% Inputs:
%    classname  - class name
%    n          - number of rows and columns
%    ni         - length of i-th dimension
%
% Outputs:
%    A - class array of specified size
%
% Example:
%     A = emptyClassArray('zonotope',5,1,2);

% Authors:       Victor Gassmann
% Written:       06-July-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% input arg checks
narginchk(1,Inf);
if nargin == 2
    if size(varargin{1},1)~=1
        throw(CORAerror('CORA:wrongValue','first','1 by x size vector'));
    end
    sz = varargin{1};
elseif nargin == 3 && any(cellfun(@(cc)numel(cc)~=1,varargin))
    throw(CORAerror('CORA:wrongValue','first/second/third/...','integer'));
else
    sz = horzcat(varargin{:});
end

inputArgsCheck({{classname,'att',{'char','string'},{'istextscalar'}};
                    {sz,'att',{'numeric'},{'integer','nonnegative'}}});

% if one value is empty, use className.empty function to construct
A = repmat(eval(classname),sz);

% ------------------------------ END OF CODE ------------------------------
