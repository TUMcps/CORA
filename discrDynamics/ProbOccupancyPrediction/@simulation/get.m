function val = get(a, propName)
% Purpose:  Get asset properties from the specified object
% Pre:      simulation object
% Post:     property value
% Built:    17.06.08, MA
% Modified: 12.10.09, MA
%           03.11.09, MA
%           04.11.09, MA

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
    error([propName,' Is not a valid asset property'])
end