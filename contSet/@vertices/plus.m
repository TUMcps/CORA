function [V] = plus(summand1,summand2)
% plus - Overloaded '+' operator for addition of a vertices object with a 
% vector
%
% Syntax:  
%    [V] = plus(summand1,summand2)
%
% Inputs:
%    summand1 - vertices object or numerical vector
%    summand2 - vertices object or numerical vector
%
% Outputs:
%    V - vertices object after addition
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      07-October-2008 
% Last update:  23-March-2015
%               24-August-2016
% Last revision:---

%------------- BEGIN CODE --------------

%Find a vertices object
%Is summand1 a vertices object?
if isa(summand1,'vertices')
    %initialize resulting vertices object
    V=summand1;
    %initialize other summand
    summand=summand2;
%Is summand2 a vertices object?    
elseif isa(summand2,'vertices')
    %initialize resulting vertices object
    V=summand2;
    %initialize other summand
    summand=summand1;  
end
    
%is summand a vector?
if isnumeric(summand)
    %Calculate sum
    for i = 1:length(V.V(1,:))
        V.V(:,i)=V.V(:,i)+summand;
    end

%something else?    
else
    V.V=[];
    disp('this operation is not implemented');
end

%------------- END OF CODE --------------