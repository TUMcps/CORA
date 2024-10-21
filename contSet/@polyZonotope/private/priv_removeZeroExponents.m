function [E,G,c] = priv_removeZeroExponents(E,G)
% priv_removeZeroExponents - add up all generators that belong to zero
%    terms
%
% Syntax:
%    [E,G,c] = priv_removeZeroExponents(E,G)
%
% Inputs:
%    E - matrix containing the exponent vectors
%    G - generator matrix
%
% Outputs:
%    E - modified exponent matrix
%    G - modified generator matrix
%    c - extracted center
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann
% Written:       12-January-2021 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% all-zero columns in exponent matrix
ind = sum(E,1)==0;
% add corresponding columns in dependent generator matrix to the center
c = sum(G(:,ind),2);
% remove corresponding columns in dependent generator matrix
G = G(:,~ind);
% remove all-zero columns in exponent matrix
E(:,ind) = [];

% ------------------------------ END OF CODE ------------------------------
