function options = rA_propInhom(obj, options)
% rA_propInhom - propagate inhomogeneous solution,
%    including reduction of zonotope order until limit given by
%    allowed bloating reached
% formula given by:
% P^R([0, t_i + Delta t_i]) = P^R([0, t_i])
%    + e^At_i * P^R_U(Delta t_i) + e^At_i-1 * P^R_u(Delta t_i-1)
%
% Syntax:  
%    options = rA_propInhom(obj,options)
%
% Inputs:
%    obj     - linearSys object
%    options - options struct
%
% Outputs:
%    options - options struct
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: -
% 
% Author:       Mark Wetzlinger
% Written:      29-Aug-2019
% Last update:  15-October-2019
% Last revision:---

%------------- BEGIN CODE --------------

% calculate linearly growing boundary for allowed bloating error
options.eSadm = options.eSmax * options.t / options.tFinal;

% additional inhomogenuity due to current time step
Rinhomadd = options.Q * options.Raux + options.P * options.lastRtrans;

if options.isInput
    % generators (+length) of added part
    options.RinhomG = [options.RinhomG, generators(Rinhomadd)];
    options.RinhomGlength = [options.RinhomGlength, vecnorm(generators(Rinhomadd),2)];

    % sort generators by length
    numGensInhom = length(options.RinhomGlength);
    [~, idxRinhomG] = mink(options.RinhomGlength,numGensInhom);
    
    binary = true; % binary search vs. start-to-end search
    
    if binary
    
    lower = 0;
    upper = numGensInhom + 1; 
    gen = floor((upper - lower) * 0.5);
    while true
        addErrorS = vecnorm(sum(abs(options.RinhomG(:,idxRinhomG(1:gen))),2),2);
        if options.eS + addErrorS > options.eSadm
            % solution < gen
            if gen == 1
                break
            end
            upper = gen;
        else
            % solution > gen
            if gen == lower
                % gen is max. number that can be over-approximated
                break
            end
            lower = gen;
        end
        % update curr
        gen = lower + floor((upper - lower) * 0.5); 
    end
    
    redIdx = idxRinhomG(1:gen);
    nonredIdx = idxRinhomG(gen+1:numGensInhom);
    options.RinhomConvInt = options.RinhomConvInt + ...
        diag(sum(abs(options.RinhomG(:,redIdx)),2));
    options.eS = norm(sum(abs(options.RinhomConvInt),2),2);
    % carry only non-reduced generators
    options.RinhomG = options.RinhomG(:,nonredIdx);
    options.RinhomGlength = options.RinhomGlength(nonredIdx);
    % add centers, take [non-reduced generators, reduced generators]
    options.Rinhom = ...
        zonotope([center(options.Rinhom)+center(Rinhomadd),...
        options.RinhomG,options.RinhomConvInt]);
    else
    
    for gen=1:numGensInhom
        addErrorS = norm(sum(abs(options.RinhomG(:,idxRinhomG(1:gen))),2),2);
        if options.eS + addErrorS > options.eSadm || gen == numGensInhom
            nonredIdx = idxRinhomG(gen:end);
            options.RinhomConvInt = options.RinhomConvInt + ...
                diag(sum(abs(options.RinhomG(:,idxRinhomG(1:gen-1))),2));
            options.eS = norm(sum(abs(options.RinhomConvInt),2),2);
            % carry only non-reduced generators
            options.RinhomG = options.RinhomG(:,nonredIdx);
            options.RinhomGlength = options.RinhomGlength(nonredIdx);
            % add centers, take [non-reduced generators, reduced generators]
            options.Rinhom = ...
                zonotope([center(options.Rinhom)+center(Rinhomadd),...
                options.RinhomG,options.RinhomConvInt]);
            break
        end
    end

    end
    
    % monitor how many generators are reduced
    options.RinhomRed(end+1,1) = gen;
end

% update Rtrans
options.lastRtrans = options.Rtrans;




end

%------------- END OF CODE --------------

