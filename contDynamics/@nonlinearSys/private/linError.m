function [error] = linError(obj,options,R)
% linError - computes the linearization error
%
% Syntax:  
%    [obj] = linError(obj,options)
%
% Inputs:
%    obj - nonlinear system object
%    options - options struct
%    R - actual reachable set
%
% Outputs:
%    obj - nonlinear system object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff
% Written:      29-October-2007 
% Last update:  22-January-2008
%               02-February-2010
%               25-July-2016 (intervalhull replaced by interval)
% Last revision: ---

%------------- BEGIN CODE --------------

%compute interval of reachable set
IH=interval(R);

%compute intervals of total reachable set
totalInt=interval(IH) + obj.linError.p.x;

%compute intervals of input
if strcmp('interval',class(options.U))
    IHinput=options.U + options.uTrans;
else
    IHinput=interval(options.U) + options.uTrans;
end
inputInt=interval(IHinput);

% %translate intervals by linearization point
% IH=IH+(-obj.linError.p.x);
IHinput=IHinput + (-obj.linError.p.u);

%obtain maximum absolute values within IH, IHinput
IHinf=abs(infimum(IH));
IHsup=abs(supremum(IH));
dx=max(IHinf,IHsup);

IHinputInf=abs(infimum(IHinput));
IHinputSup=abs(supremum(IHinput));
du=max(IHinputInf,IHinputSup);

%compute linearization error by passing the intervals to the Lagrange
%remainder mFile
if isfield(options,'lagrangeRem') && isfield(options.lagrangeRem,'method') && ...
   ~strcmp(options.lagrangeRem.method,'interval')

    % create taylor models or zoo-objects
    [objX,objU] = initRangeBoundingObjects(totalInt,inputInt,options);

    % evaluate the Lagrane remainder 
    error=obj.lagrangeRemainder(objX,objU,dx,du);
else
    error=obj.lagrangeRemainder(totalInt,inputInt,dx,du);
end
error=supremum(abs(error));

%------------- END OF CODE --------------