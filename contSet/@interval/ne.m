function res = ne( int1, int2 )
% ne( int1, int2 ) - ' ~= ' overloading
%
% Syntax:  
%    res = ne( int1, int2 )
%
% Inputs:
%    int1, int2 - intervals
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: interval
% Subfunctions: none
% MAT-files required: none
%
%
% Author:       Dmitry Grebenyuk
% Written:      06-August-2017
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    res = any(infimum(int1) ~= infimum(int2)) ||...
        any(supremum(int1) ~= supremum(int2));

end

%----------------- END -----------------