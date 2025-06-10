function res = hasTimeInterval(R)
% hasTimeInterval - checks if time interval solution is present
%
% Syntax:
%    res = hasTimeInterval(R)
%
% Inputs:
%    R - reachSet object
%
% Outputs:
%    res = true;
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachSet

% Authors:       Tobias Ladner
% Written:       30-April-2025       
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = ~isempty(R(1).timeInterval) && ...
    ~isempty(R(1).timeInterval.set);

end

% ------------------------------ END OF CODE ------------------------------
