function I = reshape(I,varargin)
% reshape - Overloads the operator 'reshape' for reshaping matrices
%
% Syntax:  
%    I = reshape(I,varargin)
%
% Inputs:
%    I - interval object
%    sz1,...,szN - integers defining the reshaping
%
% Outputs:
%    I - interval object 
%
% Example: 
%    I = interval([-1 -2; -3 -4], [1 2; 3 4]);
%    I = reshape(I,4,1);
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

%apply reshaping for infimum and supremum
I.inf = reshape(I.inf, varargin{1:end});
I.sup = reshape(I.sup, varargin{1:end});

%------------- END OF CODE --------------