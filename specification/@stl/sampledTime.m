function res = sampledTime(obj,dt,varargin)
% sampledTime - convert STL formula to sampled-time STL according to
%               Section 4.2 in [1]
%
% Syntax:
%    res = sampledTime(obj,dt)
%    res = sampledTime(obj,dt,isNNF)
%    res = sampledTime(obj,dt,isNNF,timeInt)
%    res = sampledTime(obj,dt,isNNF,timeInt,timePoint)
%
% Inputs:
%    obj - logic formula (class stl)
%    dt - time step size for sampling
%    isNNF - boolean specifiying if formula is in negation normal form
%    timeInt - matrix where the entry at index (i,j) speficies that the
%              predicate with id i at time step j is satisfied (value 1), 
%              violated (value 0), or unknown (value 2) for the time
%              interval reachable set
%    timePoint - same a timeInt but for the time point reachable set
%
% Outputs:
%    res - resulting sampled time stl formula (class stl)
%
% Example: 
%    x = stl('x',2);
%    eq = ~(x(1) < 5 | globally(x(2) < 3,interval(0.1,0.2)));
%    eq_ = sampledTime(eq,0.1)
%
% References: 
%    [1] H. Roehm et al. "STL Model Checking of Continuous and Hybrid
%        Systems", International Symposium on Automated Technology for 
%        Verification and Analysis, pp. 412-427, 2016.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Authors:       Niklas Kochdumper, Benedikt Seidl
% Written:       09-November-2022 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse input arguments
    isNNF = false; timeInt = []; timePoint = [];

    if nargin > 2 && ~isempty(varargin{1})
        isNNF = varargin{1};
    end
    if nargin > 3 && ~isempty(varargin{2})
        timeInt = varargin{2};
    end
    if nargin > 4 && ~isempty(varargin{3})
        timePoint = varargin{3};
    end

    if isempty(timePoint) && isempty(timeInt)
        pred = [];
    else
        pred.timePoint = timePoint;
        pred.timeInt = timeInt;
    end

    % convert to negation normal form as required by Lemma 4 in [1]
    if ~isNNF
        obj = negationNormalForm(obj);
    end

    % get c-divisible STL formula by scaling with the time step size
    obj = aux_scaleTime(obj,dt);

    % apply the rewriting rules in Lemma 4 and Table 1 in [1]
    obj = aux_rewrite(obj,pred);

    % eliminate true/false terms introduced during the conversion
    res = eliminateTrueFalse(obj);
    
end


% Auxiliary functions -----------------------------------------------------

function res = aux_rewrite(obj,pred)
% recursive function that applies the rewriting rules in Lemma 4 and 
% Table 1 in [1] to convert from STL to sampled-time STL

    % substitude predicates that are already satisfied by the reachable set
    % by true/false
    if ~isempty(pred) && (~obj.temporal || ...
          (ismember(obj.type,{'finally','globally'}) && ~obj.lhs.temporal))
        obj = evalPredicates(obj,pred.timeInt,pred.timePoint);
    end

    % apply rewriting rules from Lemma 4 and Table 1 in [1]
    if ~obj.temporal

        res = obj;

    elseif strcmp(obj.type,'until')

        % 1st rule in Lemma 4 in [1]
        if obj.from == 0 && obj.to == 0

            res = aux_rewrite(obj.rhs,pred);

        % 2nd rule in Table 1 in [1]
        elseif obj.from == 0

            I = interval(0,1);
            I_ = interval(0,obj.to-1);
            phi1 = obj.lhs; phi2 = obj.rhs;

            eq = phi2 | (phi1 & globally(phi1,I) & ...
                        (finally(phi2,I) | next(until(phi1,phi2,I_),1)));

            res = aux_rewrite(eq,pred);

        % 1st rule in Table 1 in [1]
        else

            I = interval(0,1);
            I_ = interval(obj.from-1,obj.to-1);
            phi1 = obj.lhs; phi2 = obj.rhs;

            eq = phi1 & globally(phi1,I) & next(until(phi1,phi2,I_),1);

            res = aux_rewrite(eq,pred);
        end

    elseif strcmp(obj.type,'release')

        % 2nd rule in Lemma 4 in [1]
        if obj.from == 0 && obj.to == 0

            res = aux_rewrite(obj.rhs,pred);

        % 4th rule in Table 1 in [1]
        elseif obj.from == 0

            I = interval(0,1);
            I_ = interval(0,obj.to-1);
            phi1 = obj.lhs; phi2 = obj.rhs;

            eq = phi2 & (phi1 | (globally(phi2,I) & ... 
                    (finally(phi1,I) | next(release(phi1,phi2,I_),1))));

            res = aux_rewrite(eq,pred);

        % 3rd rule in Table 1 in [1]
        else

            I = interval(0,1);
            I_ = interval(obj.from-1,obj.to-1);
            phi1 = obj.lhs; phi2 = obj.rhs;

            eq = phi1 | finally(phi1,I) | next(release(phi1,phi2,I_),1);

            res = aux_rewrite(eq,pred);

        end

    elseif strcmp(obj.type,'finally')

        % 3rd rule in Lemma 4 in [1]
        if strcmp(obj.lhs.type,'next')
    
            I_ = interval(obj.from,obj.to);

            res = aux_rewrite(next(finally(obj.lhs.lhs,I_),1),pred);

        % 7th rule in Table 1 in [1]
        elseif strcmp(obj.lhs.type,'&') && obj.lhs.temporal && ...
                obj.from == 0 && obj.to == 1

            I = interval(0,1);
            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs;

            eq = (globally(phi1,I) & finally(phi2,I)) | ...
                                    (finally(phi1,I) & globally(phi2,I));

            res = aux_rewrite(eq,pred);

        % 8th rule in Table 1 in [1]
        elseif strcmp(obj.lhs.type,'|') && obj.lhs.temporal && ...
                obj.from == 0 && obj.to == 1

            I = interval(0,1);
            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs;

            eq = finally(phi1,I) | finally(phi2,I);

            res = aux_rewrite(eq,pred);

        % 10th rule in Table 1 in [1]
        elseif strcmp(obj.lhs.type,'until') && obj.lhs.from == 0 && ...
                obj.lhs.to > 0 && obj.from == 0 && obj.to == 1

            I = interval(0,1);
            I_ = interval(0,obj.lhs.to-1);
            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs;

            eq = finally(phi2,I) | (globally(phi1,I) & ...
                    next(phi2 | (phi1 & globally(phi1,I) & ...
                    finally(until(phi1,phi2,I_),I)),1));

            res = aux_rewrite(eq,pred);

        % 9th rule in Table 1 in [1]
        elseif strcmp(obj.lhs.type,'until') && obj.lhs.from > 0 && ...
                obj.lhs.to > 0 && obj.from == 0 && obj.to == 1

            I = interval(0,1);
            I_ = interval(obj.lhs.from-1,obj.lhs.to-1);
            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs;

            eq = globally(phi1,I) & next(phi1 & globally(phi1,I) & ...
                    finally(until(phi1,phi2,I_),I),1);

            res = aux_rewrite(eq,pred);

        % 14th rule in Table 1 in [1]
        elseif strcmp(obj.lhs.type,'release') && obj.lhs.from == 0 && ...
                obj.lhs.to > 0 && obj.from == 0 && obj.to == 1

            I = interval(0,1);
            I_ = interval(0,obj.lhs.to-1);
            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs;

            eq = finally(phi1 & phi2,I) | (globally(phi2,I) & ...
                   next(phi2 & (phi1 | (globally(phi2,I) & ...
                   finally(release(phi1,phi2,I_),I))),1));

            res = aux_rewrite(eq,pred);

        % 13th rule in Table 1 in [1]
        elseif strcmp(obj.lhs.type,'release') && obj.lhs.from > 0 && ...
                obj.lhs.to > 0 && obj.from == 0 && obj.to == 1

            I = interval(0,1);
            I_ = interval(obj.lhs.from-1,obj.lhs.to-1);
            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs;

            eq = finally(phi1,I) | next(phi1 | finally(phi1,I) | ...
                                    finally(release(phi1,phi2,I_),I),1); 

            res = aux_rewrite(eq,pred);

        elseif obj.from ~= 0 || obj.to ~= 1

            res = aux_rewrite(until(true,obj.lhs, ...
                                interval(obj.from,obj.to)),pred);

        elseif strcmp(obj.lhs.type,'finally') && obj.lhs.from == 0 && ...
                obj.lhs.to == 1 && obj.from == 0 && obj.to == 1

            I1 = interval(obj.from,obj.to);
            I2 = interval(obj.lhs.from,obj.lhs.to);

            res = aux_rewrite(finally(until(true,obj.lhs.lhs,I2),I1),pred);

        elseif strcmp(obj.lhs.type,'globally') && obj.lhs.from == 0 && ...
                obj.lhs.to == 1 && obj.from == 0 && obj.to == 1

            I1 = interval(obj.from,obj.to);
            I2 = interval(obj.lhs.from,obj.lhs.to);

            res = aux_rewrite(finally(release(false,obj.lhs.lhs,I2),I1),pred);

        elseif obj.lhs.temporal

            I_ = interval(obj.from,obj.to);

            res = aux_rewrite(finally(aux_rewrite(obj.lhs,pred),I_),pred);

        else
            res = obj;
        end

    elseif strcmp(obj.type,'globally')

        % 4th rule in Lemma 4 in [1]
        if strcmp(obj.lhs.type,'next')
    
            I_ = interval(obj.from,obj.to);

            res = aux_rewrite(next(globally(obj.lhs.lhs,I_),1),pred);

        % 5th rule in Table 1 in [1]
        elseif strcmp(obj.lhs.type,'&') && obj.lhs.temporal && ...
               obj.from == 0 && obj.to == 1

            I = interval(0,1);
            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs;

            eq = globally(phi1,I) & globally(phi2,I);

            res = aux_rewrite(eq,pred);

        % 6th rule in Table 1 in [1]
        elseif strcmp(obj.lhs.type,'|') && obj.lhs.temporal && ...
               obj.from == 0 && obj.to == 1

            I = interval(0,1);
            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs;

            eq = globally(phi1,I) | globally(phi2,I);

            res = aux_rewrite(eq,pred);

        % 12th rule in Table 1 in [1]
        elseif strcmp(obj.lhs.type,'until') && obj.lhs.from == 0 && ...
                obj.lhs.to > 0 && obj.from == 0 && obj.to == 1

            I = interval(0,1);
            I_ = interval(0,obj.lhs.to-1);
            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs;

            eq = globally(phi2,I) | globally(phi1,I) & ...
                    next(phi2 | (phi1 & globally(phi1,I) & ...
                            globally(until(phi1,phi2,I_),I)),1);

            res = aux_rewrite(eq,pred);

        % 11th rule in Table 1 in [1]
        elseif strcmp(obj.lhs.type,'until') && obj.lhs.from > 0 && ...
                obj.lhs.to > 0 && obj.from == 0 && obj.to == 1

            I = interval(0,1);
            I_ = interval(obj.lhs.from-1,obj.lhs.to-1);
            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs;

            eq = globally(phi1,I) & next(phi1 & globally(phi1,I) & ...
                                    globally(until(phi1,phi2,I_),I),1);

            res = aux_rewrite(eq,pred);

        % 16th rule in Table 1 in [1]
        elseif strcmp(obj.lhs.type,'release') && obj.lhs.from == 0 && ...
                obj.lhs.to > 0 && obj.from == 0 && obj.to == 1

            I = interval(0,1);
            I_ = interval(0,obj.lhs.to-1);
            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs;

            eq = globally(phi2,I) & globally(phi1,I) | ...
                   next(phi2 & (phi1 | (globally(phi2,I) & ...
                   globally(release(phi1,phi2,I_),I))),1);

            res = aux_rewrite(eq,pred);
    
        % 15th rule in Table 1 in [1]
        elseif strcmp(obj.lhs.type,'release') && obj.lhs.from > 0 && ...
                obj.lhs.to > 0 && obj.from == 0 && obj.to == 1

            I = interval(0,1);
            I_ = interval(obj.lhs.from-1,obj.lhs.to-1);
            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs;

            eq = globally(phi1,I) | next(phi1 | ...
                                        globally(release(phi1,phi2,I_)));

            res = aux_rewrite(eq,pred);
        
        elseif obj.from ~= 0 || obj.to ~= 1

            res = aux_rewrite(release(false,obj.lhs, ...
                            interval(obj.from,obj.to)),pred);

        elseif strcmp(obj.lhs.type,'globally') && obj.lhs.from == 0 && ...
                obj.lhs.to == 1 && obj.from == 0 && obj.to == 1

            I1 = interval(obj.from,obj.to);
            I2 = interval(obj.lhs.from,obj.lhs.to);

            res = aux_rewrite(globally(release(false,obj.lhs.lhs,I2),I1),pred);

        elseif strcmp(obj.lhs.type,'finally') && obj.lhs.from == 0 && ...
                obj.lhs.to == 1 && obj.from == 0 && obj.to == 1

            I1 = interval(obj.from,obj.to);
            I2 = interval(obj.lhs.from,obj.lhs.to);

            res = aux_rewrite(globally(until(true,obj.lhs.lhs,I2),I1),pred);

        elseif obj.lhs.temporal

            I_ = interval(obj.from,obj.to);

            res = aux_rewrite(globally(aux_rewrite(obj.lhs,pred),I_),pred);

        else
            res = obj;
        end

    elseif strcmp(obj.type,'next')

        % shift predicate evaluations
        if ~isempty(pred)
            if ~isempty(pred.timePoint)
                pred.timePoint = pred.timePoint(:,obj.from+1:end);
            end
            if ~isempty(pred.timeInt)
                pred.timeInt = pred.timeInt(:,obj.from+1:end);
            end
        end

        % apply rewriting rules
        inner = aux_rewrite(obj.lhs,pred);

        if strcmp(inner.type,'true') || strcmp(inner.type,'false')
            res = inner;
        else
            res = next(inner,obj.from);
        end

    elseif strcmp(obj.type,'&')

        lhs = aux_rewrite(obj.lhs,pred);

        if strcmp(lhs.type,'false')
            res = stl(false);
        elseif strcmp(lhs.type,'true')
            res = aux_rewrite(obj.rhs,pred);
        else
            rhs = aux_rewrite(obj.rhs,pred);

            if strcmp(rhs.type,'false')
                res = stl(false);
            elseif strcmp(rhs.type,'true')
                res = lhs;
            else
                res = lhs & rhs;
            end
        end

    elseif strcmp(obj.type,'|')

        lhs = aux_rewrite(obj.lhs,pred);

        if strcmp(lhs.type,'true')
            res = stl(true);
        elseif strcmp(lhs.type,'false')
            res = aux_rewrite(obj.rhs,pred);
        else
            rhs = aux_rewrite(obj.rhs,pred);

            if strcmp(rhs.type,'true')
                res = stl(true);
            elseif strcmp(rhs.type,'false')
                res = lhs;
            else
                res = lhs | rhs;
            end
        end
    end
end

function res = aux_scaleTime(obj,dt)
% recursive function to scale the time by the time step size dt

    scale = 1/dt;

    if ~obj.temporal

        res = obj;

    elseif strcmp(obj.type,'until')

        res = obj;
        from = scale*res.from; to = scale*res.to;

        if abs(round(from) - from) > 1e-10 || abs(round(to) - to) > 1e-10
            res.from = ceil(from); res.to = floor(to);
        else
            res.from = round(from); res.to = round(to);
        end

        res.lhs = aux_scaleTime(res.lhs,dt);
        res.rhs = aux_scaleTime(res.rhs,dt);

    elseif strcmp(obj.type,'release')

        res = obj;
        from = scale*res.from; to = scale*res.to;

        if abs(round(from) - from) > 1e-10 || abs(round(to) - to) > 1e-10
            res.from = floor(from); res.to = ceil(to);
        else
            res.from = round(from); res.to = round(to);
        end

        res.lhs = aux_scaleTime(res.lhs,dt);
        res.rhs = aux_scaleTime(res.rhs,dt);

    elseif strcmp(obj.type,'globally')

        res = obj;
        from = scale*res.from; to = scale*res.to;

        if abs(round(from) - from) > 1e-10 || abs(round(to) - to) > 1e-10
            res.from = floor(from); res.to = ceil(to);
        else
            res.from = round(from); res.to = round(to);
        end

        res.lhs = aux_scaleTime(res.lhs,dt);

    elseif strcmp(obj.type,'finally')

        res = obj;
        from = scale*res.from; to = scale*res.to;

        if abs(round(from) - from) > 1e-10 || abs(round(to) - to) > 1e-10
            res.from = ceil(from); res.to = floor(to);
        else
            res.from = round(from); res.to = round(to);
        end

        res.lhs = aux_scaleTime(res.lhs,dt);

    elseif strcmp(obj.type,'next')

        res = obj;
        from = scale*res.from;

        if abs(round(from) - from) > 1e-10 
            inner = aux_scaleTime(res.lhs,dt);
            res = globally(inner,interval(floor(from),ceil(from)));
        else
            res.from = round(from);
            res.lhs = aux_scaleTime(res.lhs,dt);
        end

    elseif strcmp(obj.type,'~')
        
        res = ~aux_scaleTime(obj.lhs,dt);

    else

        res = obj;
        res.lhs = aux_scaleTime(res.lhs,dt);
        res.rhs = aux_scaleTime(res.rhs,dt);
    end
end

% ------------------------------ END OF CODE ------------------------------
