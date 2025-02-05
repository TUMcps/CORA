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
    case 'prob'
        val = a.result.p; 
    case 'total'
        val = a.result.pTotal;     
    case 'posProb'
        val = a.result.positionProbability; 
    case 'velProb'
        val = a.result.velocityProbability;
    case 'posProb_T'
        val = a.result.positionProbability_T; 
    case 'velProb_T'
        val = a.result.velocityProbability_T;    
    case 'inputProb'
        val = a.result.inputProbability;      
    case 'avgVel'
        val = a.result.avgVelocity;          
    case 'lcEvolProb'
        val = a.result.lcEvolProb;        
otherwise
    throw(CORAerror('CORA:specialError',[propName,' Is not a valid asset property']))
end

% ------------------------------ END OF CODE ------------------------------
