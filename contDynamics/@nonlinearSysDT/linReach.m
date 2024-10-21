function [Rtp,options,Verror] = linReach(nlnsysDT,Rinit,params,options)
% linReach - computes the reachable set after linearization
%
% Syntax:
%    [Rtp,options,Verror] = linReach(nlnsysDT,Rinit,params,options)
%
% Inputs:
%    nlnsysDT - nonlinearSysDT object
%    Rinit - initial reachable set
%    params - model parameters
%    options - options struct
%
% Outputs:
%    Rtp - resulting reachable set
%    options - options struct
%    Verror - linearization error
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Niklas Kochdumper, Mark Wetzlinger
% Written:       21-August-2012
% Last update:   29-January-2018 (NK)
%                08-April-2021 (NK, use exact plus for polyZonotopes)
%                18-June-2021 (MW, adaptive algorithm)
%                06-November-2023 (LL, add Verror as output)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % linearize nonlinear system
    [nlnsysDT,A_lin,U] = linearize(nlnsysDT,Rinit,params); 

    % translate Rinit by linearization point
    Rdelta = Rinit + (-nlnsysDT.linError.p.x);

    % compute reachable set of linearized system
    if (isa(Rdelta,'polyZonotope') || isa(Rdelta,'conPolyZono')) && ...
            (isa(U,'polyZonotope') || isa(U,'conPolyZono'))
        Rtp = exactPlus(A_lin*Rdelta, U);
    else
        Rtp = A_lin*Rdelta + U;
    end
    
    % first step of adaptive: decide tensorOrder (also compute Verror)
    Verror = 0;
    if strcmp(options.alg,'lin-adaptive') && options.i == 1
        [options,Verror] = aux_tuneTensorOrder(nlnsysDT,Rdelta,params,options,[],[]);
        
    else
        % obtain abstraction error
        if options.tensorOrder == 2
            Verror = linError_mixed_noInt(nlnsysDT, Rdelta, params, options);
        elseif options.tensorOrder == 3
            Verror = linError_thirdOrder(nlnsysDT, Rdelta, params, options);
        end
    end

    % add set of abstraction errors
    if isa(Rtp,'polyZonotope') || isa(Rtp,'conPolyZono')   
        Rtp = exactPlus(Rtp,Verror);
    else
        Rtp = Rtp + Verror;
    end
    
    % zonotope order reduction
    if contains(options.alg,'adaptive')
        Rtp = reduce(Rtp,'adaptive',options.redFactor);
    else
        Rtp = reduce(Rtp,options.reductionTechnique,options.zonotopeOrder);
    end
    
    
    % only polyZonotope: restructure the polynomial zonotope
    if isa(Rtp,'polyZonotope') && isfield(options,'polyZono') && ...
        ~isinf(options.polyZono.maxPolyZonoRatio)

        temp = options.polyZono; 
        ratio = approxVolumeRatio(Rtp,temp.volApproxMethod);

        if ratio > temp.maxPolyZonoRatio
            Rtp = restructure(Rtp,temp.restructureTechnique,...
                                temp.maxDepGenOrder);
        end
    end
    
    
    % only adaptive: decide abstraction order for next step
    if strcmp(options.alg,'lin-adaptive') && options.i > 1
        radVerror = rad(interval(Verror));
        if options.tensorOrder == 2 && ...
                all( 1 - radVerror ./ options.Verrorprev > 1-options.zetaK)
        	options = aux_tuneTensorOrder(nlnsysDT,Rdelta,params,options,radVerror,[]);
        elseif options.tensorOrder == 3 && ...
                all( 1 - radVerror ./ options.Verrorprev < options.zetaK-1)
            options = aux_tuneTensorOrder(nlnsysDT,Rdelta,params,options,[],radVerror);
        end
    end
    
end


% Auxiliary functions -----------------------------------------------------

function [options,Verror] = aux_tuneTensorOrder(nlnsysDT,Rdelta,params,options,radVerror_2,radVerror_3)

    % 1. compute other order (Verror_2 or Verror_3) if necessary
    % 2. compare to current Verror
    % 3. adapt tensorOrder and Verrorprev if change significant enough

    if isempty(radVerror_2)
        Verror_2 = linError_mixed_noInt(nlnsysDT, Rdelta, params, options);
        radVerror_2 = rad(interval(Verror_2));
    end
    if isempty(radVerror_3)
        Verror_3 = linError_thirdOrder(nlnsysDT, Rdelta, params, options);
        radVerror_3 = rad(interval(Verror_3));
    end

    if all( radVerror_3 ./ radVerror_2 > options.zetaK )
        options.tensorOrder = 2;
        if options.i == 1
            Verror = Verror_2;
        end
        options.Verrorprev = radVerror_2;
    else
        options.tensorOrder = 3;
        if options.i == 1
            Verror = Verror_3;
        end
        options.Verrorprev = radVerror_3;
    end

end

% ------------------------------ END OF CODE ------------------------------
