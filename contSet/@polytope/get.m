function val = get(P,propName)
% get - Retrieve object data from object
%
% Syntax:
%    val = get(P,propName)
%
% Inputs:
%    P - polytope object
%    propName - name of property
%
% Outputs:
%    val - value of property
%
% Example:
%    A = [1 2; -1 2; -2 -2; 1 -2];
%    b = ones(4,1);
%    P = polytope(A,b);
%    get(P,'A')
%    

% Authors:       Matthias Althoff
% Written:       12-February-2012
% Last update:   02-June-2012
%                29-October-2015
%                27-July-2016
%                04-April-2022 (adapted to polytope class)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
inputArgsCheck({{P,'att','polytope'}; ...
                {propName,'str',{'A','b','Ae','Be'}}});

% cases
switch propName 
    
    case 'A'
        val = P.A;
    
    case 'b'
        val = P.b;
    
    case 'Ae'
        val = P.Ae;
    
    case 'be'
        val = P.be;
        
    case 'V'
        val = P.V.val;
end

% ------------------------------ END OF CODE ------------------------------
