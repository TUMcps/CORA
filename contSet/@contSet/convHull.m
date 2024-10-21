function S_out = convHull(S1,S2,varargin)
% convHull - computes an enclosure for the convex hull of a set and
%    another set or a point
%
% Description:
%    computes the set { \lambda s_1 + (1-\lambda) s_2 | s_1,s_2 \in \mathcal{S}_1 \cup \mathcal{S}_2, \lambda \in [0,1] }
%
% Syntax:
%    S_out = convHull(S1)
%    S_out = convHull(S1,S2)
%    S_out = convHull(S1,S2,method)
%
% Inputs:
%    S1 - contSet object
%    S2 - contSet object or numeric
%    method - (optional) 'exact', 'outer', or 'inner'
%
% Outputs:
%    S_out - convex hull
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       30-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check number of input arguments
narginchk(1,3);

% nargin == 1 has a special meaning, immediately enter function
if nargin == 1
    S_out = convHull_(S1);
    return
end

% set default values
method = setDefaultValues({'exact'},varargin);

% check input arguments
inputArgsCheck({{S1,'att',{'contSet','numeric'}};
                {S2,'att',{'contSet','numeric','cell'}}; ...
                {method,'str',{'exact','outer','inner'}}});

% order input arguments according to their precendence
[S1,S2] = reorder(S1,S2);

% check dimension mismatch
equalDimCheck(S1,S2);

% call subclass method
S_out = convHull_(S1,S2,method);

% ------------------------------ END OF CODE ------------------------------
