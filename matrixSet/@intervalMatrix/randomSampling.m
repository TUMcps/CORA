function A = randomSampling(intMat,varargin)
% randomSampling - creates random samples within a matrix zonotope.
%
% Syntax:  
%    A = randomSampling(intMat,samples)
%
% Inputs:
%    intMat - interval matrix 
%    samples - number of segments
%    options - options struct
%
% Outputs:
%    A - cell array of matrices
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      22-June-2010
% Last update:  02-April-2017
%               01-November-2017 (more options added)
% Last revision:---

%------------- BEGIN CODE --------------

% set default inputs
if nargin == 1
    samples = 1;
    extreme = true;
elseif nargin == 2
    samples = varargin{1};
    extreme = true;
else
    samples = varargin{1};
    options = varargin{2};
    if isfield(options,'extremeSampling')
        extreme = options.extremeSampling;
    else
        extreme = true;
    end
end

%obtain dim, minimum and difference matrix
n = intMat.dim;
midA = center(intMat.int);
minA = infimum(intMat.int);
diffA = supremum(intMat.int)-infimum(intMat.int);

%generate sample matrices
A=cell(samples,1);
for i=1:samples
    if extreme % extreme
        A{i} = midA + 0.5*sign(2*rand(n)-1).*diffA;
    else % random
        A{i} = minA + rand(n).*diffA;
    end
end
    
%------------- END OF CODE --------------
