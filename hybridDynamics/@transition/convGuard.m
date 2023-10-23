function trans = convGuard(trans,inv,options)
% convGuard - converts the guard sets to the set representation that is
%    required for the selected guard-intersection method
%
% Syntax:
%    trans = convGuard(trans,inv,options)
%
% Inputs:
%    trans - transition object
%    inv - invariant set of the location
%    options - struct containing algorithm settings
%
% Outputs:
%    trans - modified transition object
%
% Example:
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper
% Written:       16-May-2018
% Last update:   19-June-2022 (MW, error message, use switch-case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

switch options.guardIntersect

    case 'polytope'

        if ~isa(trans.guard,'polytope')
            % convert to polytope
            trans.guard = polytope(trans.guard);
        end
        
        % intersect with invariant set
        trans.guard = and_(trans.guard,inv,'exact');

    case {'conZonotope','conZonotopeFast'}
    
        if isa(trans.guard,'polytope')
            % convert to conZonotope
            trans.guard = conZonotope(trans.guard);
        end      

    otherwise

        % throw error
        throw(CORAerror('CORA:wrongFieldValue','options.guardIntersect',...
            {'zonoGirard','hyperplaneMap','pancake','nondetGuard'}));

end

% ------------------------------ END OF CODE ------------------------------
