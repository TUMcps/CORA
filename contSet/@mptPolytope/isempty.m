function res = isempty(obj)
% is_empty - returns 1 if a pplPolytope is empty and 0 otherwise
%
% Syntax:  
%    res = is_empty(obj)
%
% Inputs:
%    obj - pplPolytope object
%
% Outputs:
%   res - result in {0,1}
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      03-February-2011
% Last update:  16-July-2015
%               20-August-2015
% Last revision:---

%------------- BEGIN CODE --------------


try %MPT 3
    res = all(isEmptySet(obj.P));
catch %MPT2
    
    %get properties
    Array = get(obj.P,'Array');
    H = get(obj.P,'H');
    K = get(obj.P,'K');
    RCheb = get(obj.P,'RCheb');

    %check for emtiness; see display function of mpt polytope code
    if isempty(Array) & H==1 & K==-Inf & RCheb==-Inf
        res=1;
    else
        res=0;
    end
end

%------------- END OF CODE --------------