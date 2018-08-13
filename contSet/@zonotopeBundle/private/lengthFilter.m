function [Gred]=lengthFilter(G,rem)
% lengthFilter - filters out short generators
%
% Syntax:  
%    [Gred]=lengthFilter(G,rem)
%
% Inputs:
%    G - matrix of generators
%    rem - number of remaining generators
%
% Outputs:
%    Gred - reduced matrix of generators
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

%prefilter generators
for i=1:length(G(1,:))
    h(i)=norm(G(:,i));
end
[value,index]=sort(h);
Gred=G(:,index((end-rem+1):end));

%------------- END OF CODE --------------
