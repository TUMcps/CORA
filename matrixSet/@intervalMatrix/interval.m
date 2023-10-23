function I = interval(intMat)
% interval - converts an interval matrix object to an interval object
%
% Syntax:
%    I = interval(intMat)
%
% Inputs:
%    intMat - intervalMatrix object
%
% Outputs:
%    I - interval object
%
% Example:
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       21-June-2010
% Last update:   25-July-2016 (intervalhull replaced by interval)
% Last revision: 18-June-2023 (MW, harmonize with other conversion methods)

% ------------------------------ BEGIN CODE -------------------------------

I = intMat.int;

% ------------------------------ END OF CODE ------------------------------
