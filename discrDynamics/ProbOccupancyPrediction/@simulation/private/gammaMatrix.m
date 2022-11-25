function [Gamma]=gammaMatrix(markovChainSpec,simOptions)
% gammaMatrices - Generates the Gamma matrix
%
% Syntax:  
%    [Gamma]=gammaMatrices(markovChainSpec)
%
% Inputs:
%    markovChainSpec - Markov-Chain specifications
%
% Outputs:
%    Gamma - cell array of Gamma matrices
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 03-July-2008 
% Last update: 22-August-2008
% Last revision: ---

%------------- BEGIN CODE --------------

%obtain number of segments
nrOfInputs=markovChainSpec.nrOfInputs;

a=simOptions.gamma;
%a=0.001;
%a=1000;

%build Gamma matrix
for i=1:nrOfInputs
    for j=1:nrOfInputs
        Gamma(i,j)=1/((i-j)^2+a);
    end
end
Gamma=normalizeMatrix(Gamma);

%------------- END OF CODE --------------