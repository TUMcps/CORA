function options = check_restructureTechnique(options,obj)
% check_restructureTechnique - checks if options.restructureTechnique
%  1) exists
%  2) takes an allowed value
%
% Syntax:
%    options = check_restructureTechnique(options,obj)
%
% Inputs:
%    options - options for object
%    obj     - system object
%
% Outputs:
%    options - updated options for object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      02-January-2020
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

    strct = 'options.polyZono';
    option = 'restructureTechnique';

   if isfield(options.polyZono,option)

       % parse string for the method
       method = options.polyZono.restructureTechnique;
       if startsWith(method,'zonotope')
            redMeth = method(9:end);
       elseif startsWith(method,'reduceFull')
            redMeth = method(11:end);
       elseif startsWith(method,'reduce')
            redMeth = method(7:end);
       else
            error(printOptionOutOfRange(obj,option,strct));
       end

       redMeth(1) = lower(redMeth(1));

       % compare to zonotope reduction techniques
       validRedTech = {'girard','combastel','pca','methA', ....
                       'methB','methC','methD','methE', ...
                       'methF','redistribute','cluster', ...
                       'scott','constOpt'};

       if ~ismember(redMeth,validRedTech)
           error(printOptionOutOfRange(obj,option,strct));
       end

   else
       options.polyZono.restructureTechnique = 'reduceGirard';
   end
end

%------------- END OF CODE --------------