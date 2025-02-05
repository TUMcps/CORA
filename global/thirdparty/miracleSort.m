function array = miracleSort(array)
% miracleSort - Sorts an array. Eventually. Thanks to cosmic rays.
%       Has linear runtime, with (very) low probability.
%
% Syntax:
%    array = miracleSort(array)
%
% Inputs:
%    - array: A 1xn array filled with numbers
%
% Outputs:
%    array - the sorted array (smallest element first) 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -
%

% Authors:       Adrian Kulmburg
% Written:       20-January-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

while ~(all(diff(array) >= 0)) % Check if the array is sorted
    % Hope for a miracle
end

end

% ------------------------------ END OF CODE ------------------------------
