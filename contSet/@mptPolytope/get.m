function val = get(obj, propName)
% get - Retrieve object data from obj
%
% Syntax:  
%    val = get(obj, propName)
%
% Inputs:
%    obj      - mptPolytope object
%    propName - name of property
%
% Outputs:
%    val - value of property
%
% Properties:
%    intervals - intervals of interval object

% Author:       Matthias Althoff
% Written:      12-February-2012
% Last update: 	02-June-2012
%               29-Octber-2015
%               27-July-2016
% Last revision:---

%------------- BEGIN CODE --------------

switch propName 
    case 'P'
        val = obj.P;
    case {'H','A'}
        try %MPT V3
            val = obj.P.A;
        catch %MPT V2
            [val,K] = double(obj.P);
        end
    case {'K','b'}
        try %MPT V3
            val = obj.P.b;
        catch %MPT V2
            [H,val] = double(obj.P);
        end
    case 'Ae'
        val = obj.P.Ae; %MPT V3 ONLY
    case 'be'
        val = obj.P.be; %MPT V3 ONLY
    otherwise
        error([propName,' is not a valid asset property'])
end

%------------- END OF CODE --------------