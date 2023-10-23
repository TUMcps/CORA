function S = initEmptySet(type)
% initEmptySet - instantiates an empty set of a contSet class
%
% Syntax:
%    S = initEmptySet(type)
%
% Inputs:
%    type - contSet class name
%
% Outputs:
%    S - instantiated set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       24-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% all supported classes
admissibleTypes = {
    'capsule','conPolyZono','conHyperplane','conZonotope','ellipsoid',...
    'halfspace','interval','levelSet','polytope','polyZonotope',...
    'probZonotope','zonoBundle','zonotope',... % contSet classes
    'emptySet','fullspace'}; % future contSet classes

% check input argument
inputArgsCheck({{type,'str',admissibleTypes}});


switch type

    case 'capsule'
        S = capsule();

    case 'conPolyZono'
        S = conPolyZono();

    case 'conHyperplane'
        S = conHyperplane();

    case 'conZonotope'
        S = conZonotope();

    case 'ellipsoid'
        S = ellipsoid();

    case 'halfspace'
        S = halfspace();

    case 'interval'
        S = interval();

    case 'levelSet'
        S = levelSet();

    case 'polytope'
        S = polytope();

    case 'polyZonotope'
        S = polyZonotope();

    case 'probZonotope'
        S = probZonotope();

    case 'zonoBundle'
        S = zonoBundle();

    case 'zonotope'
        S = zonotope();

    case 'emptySet'
        S = emptySet();

    case 'fullspace'
        S = fullspace();

end

% ------------------------------ END OF CODE ------------------------------
