function options = check_polyZono(options, obj)
% check_polyZono - checks if options.polyZono
%  1) takes an allowed value
%
% Syntax:
%    check_polyZono(options, obj)
%
% Inputs:
%    options - options for object
%    obj     - system object
%
% Outputs:
%    options - options for object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   -

% Author:       Niklas Kochdumper
% Written:      22-May-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% check all polyZono settings
if isfield(options, 'polyZono')
    
   if ~isa(options.R0,'polyZonotope') || ~strcmp(options.alg,'poly')
       warning('options.polyZono is redundant!');
   end
    
   % check setting 'maxPolyZonoRatio'
   options = check_maxPolyZonoRatio(options, obj);
   
   % check setting 'maxDepGenOrder'
   options = check_maxDepGenOrder(options, obj);
   
   % check setting 'restructureTechnique'
   options = check_restructureTechnique(options, obj);
    
end