function obj = convGuard(obj,inv,options)
% convGuard - convert the guard sets to the set representation that is
%             required for the selected guard-intersection method
%
% Syntax:  
%    obj = convGuard(obj,inv,options)
%
% Inputs:
%    obj - transition object
%    inv - invariant set of the location
%    options - struct containing algorithm settings
%
% Outputs:
%    obj - modified transition object
%
% Example: 
%
% Other m-files required: not specified
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      16-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

     if strcmp(options.guardIntersect,'polytope')
         
        if ~isa(obj.guard,'mptPolytope')
            
           % convert to mptPolytope
           obj.guard = polytope(obj.guard);          
          
        end
        
        % intersect with invariant set
        obj.guard = obj.guard & inv;
          
        
     elseif ~ismember(options.guardIntersect,{'zonoGirard','hyperplaneMap'})
        error('Wrong value for setting options.guardIntersect!');  
     end
end

%------------- END OF CODE --------------