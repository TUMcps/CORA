function I_out = enclose(I1,I2)
% enclose - encloses an interval and its affine transformation (calling
%           convHull_ as the operation is equivalent for intervals)
%
% Syntax:
%    I_out = enclose(I)
%
% Inputs:
%    I1 - interval object
%    I2 - interval object
%
% Outputs:
%    I2 - interval object
%
% Example: 
%    I1 = interval([1;-1], [2; 1]);
%    I2 = rand(2,2)*I1 + rand(2,1);
%    I_out = enclose(I1,I2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, zonotope/enclose

% Authors:       Maximilian Perschl
% Written:       18-November-2024 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% compute result via convHull
I_out = convHull(I1,I2);

% ------------------------------ END OF CODE ------------------------------
