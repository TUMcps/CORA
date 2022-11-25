function obj = reshape(varargin)
% reshape - Overloads the operator 'reshape' for reshaping matrices
%
% Syntax:  
%    obj = reshape(varargin)
%
% Inputs:
%    obj - interval object
%    sz1,...,szN - integers defining the reshaping
%
% Outputs:
%    obj - interval object 
%
% Example: 
%    a=interval([-1 -2; -3 -4], [1 2; 3 4]);
%    b=reshape(a, 4, 1);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      05-August-2015 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

obj = varargin{1};

%apply reshaping for infimum and supremum
obj.inf = reshape(obj.inf, varargin{2:end});
obj.sup = reshape(obj.sup, varargin{2:end});

%------------- END OF CODE --------------