function [ExpMat,Gnew,c] = removeZeroExponents(ExpMat,G)
% removeZeroExponents - add up all generators that belong to zero terms
%
% Syntax:  
%    [ExpMatNew,Gnew,c] = removeZeroExponents(ExpMat,G)
%
% Inputs:
%    ExpMat - matrix containing the exponent vectors
%    G - generator matrix
%
% Outputs:
%    ExpMat - modified exponent matrix
%    G - modified generator matrix
%    c - extracted center
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:        Victor Gassmann
% Written:       12-January-2021 
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

% all-zero columns in exponent matrix
ind = sum(ExpMat,1)==0;
% add corresponding columns in dependent generator matrix to the center
c = sum(G(:,ind),2);
% remove corresponding columns in dependent generator matrix
Gnew = G(:,~ind);
% remove all-zero columns in exponent matrix
ExpMat(:,ind) = [];

%------------- END OF CODE --------------