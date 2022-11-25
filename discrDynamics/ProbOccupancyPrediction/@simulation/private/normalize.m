function [pNorm]=normalize(p1,p2)
% normalize - normalizes a probability vector p1, such that its sum equals
% the sum of p2
%
% Syntax:  
%    [pNorm]=normalize(p1,p2)
%
% Inputs:
%    p1 - probability vector
%    p2 - probability vector
%
% Outputs:
%    pNorm - normalized probability vector
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 19-November-2007 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------
        
%normalize
s1=sum(p1(2:end)); %remove probabilities from outside and add them to others
s2=sum(p2);
s2/s1
if s1~=0
    pNorm=s2/s1*p1;
    pNorm(1,1)=0;
else
    pNorm=p1;
end


%------------- END OF CODE --------------