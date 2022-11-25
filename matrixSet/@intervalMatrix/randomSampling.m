function A = randomSampling(intA,varargin)
% randomSampling - creates random samples within a matrix zonotope.
%
% Syntax:  
%    A = randomSampling(intA,samples)
%
% Inputs:
%    intA - interval matrix 
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
    extreme = 1;
elseif nargin == 2
    samples = varargin{1};
    extreme = 1;
else
    samples = varargin{1};
    options = varargin{2};
    if isfield(options,'extremeSampling')
        extreme = options.extremeSampling;
    else
        extreme = 1;
    end
end


%obtain dim, minimum and difference matrix
dim = intA.dim;
midA = center(intA.int);
minA = infimum(intA.int);
diffA = supremum(intA.int)-infimum(intA.int);

%generate sample matrices
A=cell(samples,1);
for i=1:samples
    % extreme
    if extreme
        A{i}=midA + 0.5*sign(2*rand(dim)-1).*diffA;
    % random
    else
        A{i}=minA + rand(dim).*diffA;
    end
end

    
%------------- END OF CODE --------------
