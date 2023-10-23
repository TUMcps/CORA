function varargout = size(I, varargin)
% size - Overloads the operator that returns the size of the object
%
% Syntax:
%    varargout = size(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    varargout - varying outputs, see 'doc size'
%
% Example: 
%    I = interval([-1 1], [1 2]);
%    size(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       26-January-2016 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%return size of infimum
if nargin > 1
    [varargout{1:nargout}] = size(I.inf,varargin{1});
else
    [varargout{1:nargout}] = size(I.inf);
end

% ------------------------------ END OF CODE ------------------------------
