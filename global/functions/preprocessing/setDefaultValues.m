function varargout = setDefaultValues(defaultValues,varargin)
% setDefaultValues - set default values for input arguments; if some value
%    is given, then this value is kept (check for admissible values is done
%    using inputArgsCheck)
%
% Syntax:  
%    varargout = setDefaultValues(defaultValues,varargin)
%
% Inputs:
%    defaultValues - default values in a cell-array: {dV1,dV2,...,dVend}
%
% Outputs:
%    varargout - values for input arguments (1xn cell-array), where the
%                the non-given values are set to default values
%
% Example:
%    defValues = {1,'standard','lower',true};
%    givenValues = {10,'gaussian'};
%    [N,type,mode,flag] = setDefaultValues(defValues,givenValues{:});
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: inputArgsCheck

% Author:       Mingrui Wang, Mark Wetzlinger
% Written:      30-May-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% number of given input arguments
n_given = length(varargin);

% number of default values
n_default = max(size(defaultValues));

% assign default values if corresponding values are not provided
varargout = [varargin(1:n_given), defaultValues(n_given+1:n_default)];

%------------- END OF CODE --------------