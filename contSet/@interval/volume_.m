function vol = volume_(I,varargin)
% volume_ - Computes volume of an interval
%
% Syntax:
%    vol = volume_(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    vol - volume
%
% Example: 
%    I = interval([1; -1], [2; 1]);
%    vol = volume(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/volume

% Authors:       Matthias Althoff
% Written:       24-July-2016 
% Last update:   18-August-2022 (MW, include standardized preprocessing)
%                04-December-2023 (MW, fix degenerate-unbounded case)
% Last revision: 27-March-2023 (MW, rename volume_)

% ------------------------------ BEGIN CODE -------------------------------

% compute half of the diameter
r = rad(I);

if isempty(r) || ~all(r)
    % empty or degenerate
    vol = 0;
else
    % simple volume formula
    vol = prod(2*r);
end

% ------------------------------ END OF CODE ------------------------------
