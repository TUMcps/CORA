function obj = convGuard(obj,options)
% convGuard - convert all guard sets to the set representation that is
%             required for the selected guard-intersection method
%
% Syntax:  
%    obj = convGuard(obj,inv,options)
%
% Inputs:
%    obj - location object
%    options - struct containing algorithm settings
%
% Outputs:
%    obj - modified location object
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

   % convert invariant set to mptPolytope
   if ~options.isHyperplaneMap
       if (isfield(options,'intersectInvariant') && options.intersectInvariant==1) || ...
           strcmp(options.guardIntersect,'polytope')

            try
               obj.invariant = polytope(obj.invariant,options); 
            catch
               obj.invariant = polytope(obj.invariant);
            end
       end

       % convert all guard sets to the required set representation
       for i = 1:length(obj.transition)
           obj.transition{i} = convGuard(obj.transition{i},obj.invariant,options);
       end
   end
end

%------------- END OF CODE --------------