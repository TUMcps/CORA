function [G]=nonzeroFilter(G)
% nonZeroFilter - filters out generators of length 0
%
% Syntax:  
%    [G]=nonzeroFilter(G)
%
% Inputs:
%    G - matrix of generators
%
% Outputs:
%    G - reduced matrix of generators
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author: Matthias Althoff
% Written: 12-September-2008
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%Delete zero-generators
i=1;
while i<=length(G(1,:))
    if G(:,i)==0*G(:,i)
        G(:,i)=[];
    else
        i=i+1;
    end
end

%------------- END OF CODE --------------
