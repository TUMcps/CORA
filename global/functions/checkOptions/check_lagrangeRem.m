function check_lagrangeRem(options, obj)
% check_lagrangeRem - checks if options.lagrangeRem
%  1) takes an allowed value
%
% Syntax:
%    check_lagrangeRem(options, obj)
%
% Inputs:
%    options - options for object
%    obj     - system object
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   -

% Author:       Mark Wetzlinger
% Written:      04-Mar-2019
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

strct = 'options';

% check all lagrange remainder settings
if isfield(options, 'lagrangeRem')
    
    temp = options.lagrangeRem;
    
    % check setting 'method'
    if isfield(temp,'method') && ~ismember(temp.method,{'interval','taylorModel','zoo'})
        error(printOptionOutOfRange(obj,'lagrangeRem.method',strct));
    end
    
    % check setting 'zooMethods'
    if isfield(temp,'zooMethods') 
       if ~isfield(temp,'method') || ~strcmp(temp.method,'zoo')
           warning('The following params/options have been set, but are redundant:\n  lagrangeRem.zooMethods');
       else
           wrong = 0;
           values = {'interval','affine(int)','affine(bnb)', ...
                     'affine(bnbAdv)','affine(linQuad)','taylm(int)', ...
                     'taylm(bnb)','taylm(bnbAdv)','taylm(linQuad)'};
           
           if ~iscell(temp.zooMethods)
              wrong = 1; 
           else
               for i = 1:length(temp.zooMethods)
                  if ~ischar(temp.zooMethods{i}) || ~ismemeber(temp.zooMethods{i},values)
                     wrong = 1;
                     break;
                  end
               end
           end
           
           if wrong
               error(printOptionOutOfRange(obj,'lagrangeRem.zooMethods',strct));
           end
       end 
    end
    
    % check setting 'optMethod'
    if isfield(temp,'optMethod')
        if ~isfield(temp,'method') || ~strcmp(temp.method,'taylorModel')
           warning('The following params/options have been set, but are redundant:\n  lagrangeRem.optMethod');
        else
           values = {'int','bnb','bnbAdv','linQuad'};  
            
           if ~ischar(temp.optMethod) || ~ismember(temp.optMethod,values)
               error(printOptionOutOfRange(obj,'lagrangeRem.optMethod',strct));
           end
        end
    end
    
    % check setting 'maxOrder'
    if isfield(temp,'maxOrder')
        if ~isfield(temp,'method') || ~ismember(temp.method,{'taylorModel','zoo'})
           warning('The following params/options have been set, but are redundant:\n  lagrangeRem.maxOrder');
        else
           if ~isscalar(temp.maxOrder) || mod(temp.maxOrder,1) ~= 0 || temp.maxOrder <= 0
               error(printOptionOutOfRange(obj,'lagrangeRem.maxOrder',strct));
           end
        end
    end
    
    % check setting 'tolerance'
    if isfield(temp,'tolerance')
        if ~isfield(temp,'method') || ~ismember(temp.method,{'taylorModel','zoo'})
           warning('The following params/options have been set, but are redundant:\n  lagrangeRem.tolerance');
        else
           if ~isscalar(temp.tolerance) || temp.tolerance <= 0
               error(printOptionOutOfRange(obj,'lagrangeRem.tolerance',strct));
           end
        end
    end
    
    % check setting 'eps'
    if isfield(temp,'eps')
        if ~isfield(temp,'method') || ~ismember(temp.method,{'taylorModel','zoo'})
           warning('The following params/options have been set, but are redundant:\n  lagrangeRem.eps');
        else
           if ~isscalar(temp.eps) || temp.eps <= 0
               error(printOptionOutOfRange(obj,'lagrangeRem.eps',strct));
           end
        end
    end
    
    % check setting 'simplify'
    check_simplify(temp, obj);
    
    % check setting 'replacements'
    check_replacements(temp, obj);
    
    % check setting 'tensorParallel'
    check_tensorParallel(temp, obj);
end

end

%------------- END OF CODE --------------

