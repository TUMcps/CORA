function Rtp = linReach(obj,Rinit,options)
% linReach - computes the reachable set after linearization
%
% Syntax:  
%    [Rtp] = linReach(obj,Rinit,options)
%
% Inputs:
%    obj - nonlinearSysDT system object
%    Rinit - initial reachable set
%    options - options struct
%
% Outputs:
%    Rtp - resulting reachable set
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      21-August-2012
% Last update:  29-January-2018 (NK)
%               08-April-2021 (NK, use exact plus for polyZonotopes)
% Last revision:---

%------------- BEGIN CODE --------------

    % linearize nonlinear system
    [obj,A_lin,U] = linearize(obj,Rinit,options); 

    % translate Rinit by linearization point
    Rdelta = Rinit + (-obj.linError.p.x);

    % compute reachable set of linearized system
    Rtp = A_lin*Rdelta + U;

    % obtain abstraction error
    if options.tensorOrder == 2
        Verror = linError_mixed_noInt(obj, options, Rdelta);
    elseif options.tensorOrder == 3
        Verror = linError_thirdOrder(obj, options, Rdelta);    
    end

    % add set of abstraction errors
    if options.tensorOrder == 3 && ...
       (isa(Rtp,'polyZonotope') || isa(Rtp,'conPolyZono'))
   
        Rtp = exactPlus(Rtp,Verror);
    else
        Rtp = Rtp + Verror;
    end
    
    % order reduction
    Rtp = reduce(Rtp,options.reductionTechnique,options.zonotopeOrder);
    
    % potentially restructure the polynomial zonotope
    if isa(Rtp,'polyZonotope') && isfield(options,'polyZono') && ...
       ~isinf(options.polyZono.maxPolyZonoRatio)

        temp = options.polyZono; 
        ratio = approxVolumeRatio(Rtp,temp.volApproxMethod);

        if ratio > temp.maxPolyZonoRatio
           Rtp = restructure(Rtp,temp.restructureTechnique, ...
                             temp.maxDepGenOrder);
        end
    end
end

%------------- END OF CODE --------------