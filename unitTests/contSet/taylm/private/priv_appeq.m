function res = priv_appeq(int1, int2, eps)
% priv_appeq( int1, int2, eps ) - % approximaly equal with epsilon tolerance
%
% Syntax:
%    res = priv_appeq( int1, int2, eps )
%
% Inputs:
%    int1, int2 - intervals, or arrays of numbers
%
% Outputs:
%    res - boolean 
%
% Other m-files required: interval
% Subfunctions: none
% MAT-files required: none

% Authors:       Dmitry Grebenyuk
% Written:       07-August-2017
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isa(int1, 'interval') && isa(int2, 'interval')
    res = all( abs(infimum(int1) - infimum(int2)) <= eps )  &&...
        all( abs(supremum(int1) - supremum(int2)) <= eps );
elseif isa(int1, 'double') && isa(int2, 'double')
    res = all( abs(int1 - int2) <= eps );
end  
end

% ------------------------------ END OF CODE ------------------------------
