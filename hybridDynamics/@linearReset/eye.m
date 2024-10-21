function linReset = eye(n,varargin)
% eye - instantiates an identity reset function of dimension n
%
% Syntax:
%    linReset = linearReset.eye(n)
%    linReset = linearReset.eye(n,m)
%
% Inputs:
%    n - pre/post-state dimension
%    m - input dimension
%
% Outputs:
%    linReset - linearReset object
%
% Example: 
%    n = 3; m = 2;
%    linReset1 = linearReset.eye(n);
%    linReset2 = linearReset.eye(n,m);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       15-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check number of input arguments and set default number of inputs
narginchk(1,2);
m = setDefaultValues({1},varargin);

% check input arguments
inputArgsCheck({{n,'att','numeric',{'scalar','integer','nonnegative'}};
                {m,'att','numeric',{'scalar','integer','nonnegative'}}});

% instantiate reset function
if n == 0
    linReset = linearReset();
else
    linReset = linearReset(eye(n),zeros(n,m),zeros(n,1));
end

% ------------------------------ END OF CODE ------------------------------
