function X = sampleBox(Z,N)
% sampleBox - computes N samples uniformly distributed inside a
%    parallelotope
%
% Syntax:
%    X = sampleBox(Z,N)
%
% Inputs:
%    Z - zonotope object
%    N - number of samples
%
% Outputs:
%    X - samples (each sample a column vector)
%
% Example: 
%    Z = zonotope([0;0],[1 1; 0 1]);
%    X = sampleBox(Z,1000);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann
% Written:       16-October-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

G = generators(Z);
[n,nrGens] = size(G);
if n<nrGens || rank(G)<n
    throw(CORAerror('CORA:wrongValue','first',...
        'The zonotope has to be a parallelotope'));
end
options.method = 'achr';

%efficient since Z is box
HS = halfspace(Z);
A = HS.halfspace.H;
b = HS.halfspace.K;
X = cprnd(N,A,b,options);
X = X';

% ------------------------------ END OF CODE ------------------------------
