function res = abs(I)
% abs - returns the absolute value of an interval
%
% Syntax:
%    res = abs(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    res - interval object
%
% Example:
%    I = interval([-2;-1],[3;4]);
%    res = abs(I);
%
%    figure; hold on;
%    plot(I,[1,2],'b');
%    plot(res,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       26-June-2015
% Last update:   14-February-2015
%                12-October-2015
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init resulting interval object
res = I;

% separate computation for scalar and matrix case for efficiency

% scalar case
if isscalar(I)
    if I.sup < 0
        res.inf = abs(I.sup);
        res.sup = abs(I.inf);
    elseif I.inf > 0
        res = I;
    else
        res.inf = 0;
        res.sup = max(abs(I.inf), abs(I.sup));
    end
    
% matrix case
else
    
    % find negative indices (if infimum is greater than zero, the absolute value has no effect)
    ind = I.inf < 0 & I.sup > 0;  % For [-a, +b] case
    ind1 = I.inf < 0 & I.sup <= 0; % For [-a, -b] case
    
    res.sup(ind) = max(abs(I.inf(ind)), abs(I.sup(ind))); % order of computation matters
    res.inf(ind) = abs(0*I.inf(ind));
    
    res.sup(ind1) = abs(I.inf(ind1));
    res.inf(ind1) = abs(I.sup(ind1));
end

% ------------------------------ END OF CODE ------------------------------
