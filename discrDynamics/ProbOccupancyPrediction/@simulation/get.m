function val = get(a, propName)
% get - short description of the function
% Purpose:  Get asset properties from the specified object
% Pre:      simulation object
% Post:     property value
%
% Syntax:
%    val = get(a, propName)
%
% Inputs:
%    ???
%
% Outputs:
%    ???
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Matthias Althoff
% Written:       17-June-2008
% Last update:   12-October-2009
%                03-November-2009
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

switch propName
    % probability vector
    case 'prob'
        val = a.result.p; 
    % total probability vector
    case 'total'
        val = a.result.pTotal;    
    % probability vector of position
    case 'posProb'
        val = a.result.positionProbability; 
    % probability vector of velocity
    case 'velProb'
        val = a.result.velocityProbability;
    % probability vector of position for a specific point in time
    case 'posProb_T'
        val = a.result.positionProbability_T; 
    % probability vector of velocity for a specific point in time
    case 'velProb_T'
        val = a.result.velocityProbability_T;    
    % probability vector of inputs
    case 'inputProb'
        val = a.result.inputProbability;  
    % average velocity
    case 'avgVel'
        val = a.result.avgVelocity;      
    % lane change probability
    case 'lcEvolProb'
        val = a.result.lcEvolProb;        
otherwise
    throw(CORAerror('CORA:specialError',[propName,' Is not a valid asset property']))
end

% ------------------------------ END OF CODE ------------------------------
