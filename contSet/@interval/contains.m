function res = contains(I,S,varargin)
% contains - determines if an interval contains a set or a point
%
% Syntax:  
%    res = contains(I,S)
%
% Inputs:
%    I - interval object
%    S - contSet object or single point
%    type - 'exact' or 'approx'
%    tol - tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    I1 = interval([-1;-2],[2;3]);
%    I2 = interval([0;0],[1;1]);
%
%    contains(I1,I2)
%    contains(I1,[1;2])
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/contains

% Author:       Niklas Kochdumper
% Written:      16-May-2018
% Last update:  02-Sep-2019
%               15-November-2022 (MW, return logical array for points)
%               25-November-2022 (MW, rename 'contains')
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing
[res,vars] = pre_contains('interval',I,S,varargin{:});

% check premature exit
if res
    % if result has been found, it is stored in the first entry of var
    res = vars{1}; return
else
    I = vars{1}; S = vars{2}; tol = vars{4};
    % type is always exact
end


% init result
res = false;

% point in interval containment
if isnumeric(S)
    
    res = all( (I.inf < S | withinTol(I.inf,S)) ...
            & (I.sup > S | withinTol(I.sup,S)) , 1);

% interval in interval containment
elseif isa(S,'interval')

    if all(I.sup >= S.sup) && all(I.inf <= S.inf)
        res = true;
    end

% other set in interval containment
else

    res = contains(mptPolytope(I),S);

end

%------------- END OF CODE --------------