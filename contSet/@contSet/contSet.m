classdef contSet
% contSet - abstract superclass for continuous sets
%
% Syntax:  
%    S = contSet()
%    S = contSet(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    obj - generated contSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Matthias Althoff, Mark Wetzlinger
% Written:       02-May-2007 
% Last update:   04-May-2020 (MW, transition to classdef)
%                01-June-2022 (MW, add CORAerror)
%                22-March-2023 (MW, remove deprecated property dimension)
% Last revision: ---

%------------- BEGIN CODE --------------

properties (SetAccess = protected, GetAccess = public)
    % no properties
end

methods

    function obj = contSet(varargin)
        
        % (default constructor)
        if nargin == 0
            % no properties
            
        % If 1 argument is passed
        elseif nargin == 1
            % (copy constructor)
            if isa(varargin{1},'contSet')
                obj = varargin{1};
            end

        % Else if not enough or too many inputs are passed    
        else
            throw(CORAerror('CORA:tooManyInputArgs',1));
        end
    end

    
    % res = norm(S,varargin)
end

methods(Access = {?contSet, ?contDynamics})
%     function res = isIntersecting_(S1,S2,type,varargin)
%         throw(CORAerror("CORA:noops",S1,S2));
%     end
end

end

%------------- END OF CODE --------------