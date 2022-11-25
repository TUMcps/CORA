function A = randomSampling(matZ,varargin)
% randomSampling - creates random samples within a matrix zonotope
%
% Syntax:  
%    A = randomSampling(matZ,samples)
%
% Inputs:
%    matZ - matrix zonotope
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
% Written:      23-June-2010
% Last update:  02-April-2017
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

% generate sample matrices
A=cell(samples,1);
for i=1:samples
    %initialize random matrix
    A{i}=matZ.center;
    %add generator matrices
    for iGen=1:matZ.gens
        if extreme
            A{i}=A{i} + sign(2*rand(1)-1).*matZ.generator{iGen};
        else
            A{i}=A{i} + (2*rand(1)-1).*matZ.generator{iGen};
        end
    end
end

%------------- END OF CODE --------------
%------------- END OF CODE --------------
