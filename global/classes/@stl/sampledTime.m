function res = sampledTime(obj,dt)
% sampledTime - convert STL formula to sampled-time STL according to
%               Section 4.2 in [1]
%
% Syntax:  
%    res = sampledTime(obj,dt)
%
% Inputs:
%    obj - logic formula (class stl)
%    dt - time step size for sampling
%
% Outputs:
%    res - resulting sampled time stl formula (class stl)
%
% Example: 
%    x = stl('x',2)
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

% Author:       Niklas Kochdumper
% Written:      9-November-2022 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % convert to negation normal form as required by Lemma 4 in [1]
    obj = negationNormalForm(obj);

    % get c-divisible STL formula by scaling with the time step size
    obj = scaleTime(obj,dt);

    % apply the rewriting rules in Lemma 4 and Table 1 in [1]
    obj = rewrite(obj);

    % eliminate true/false terms introduced during the conversion
    res = eliminateTrueFalse(obj);
    
end

% Auxiliary Functions -----------------------------------------------------

function res = rewrite(obj)
% recursive function that applies the rewriting rules in Lemma 4 and 
% Table 1 in [1] to convert from STL to sampled-time STL

    I = interval(0,1);

    if ~obj.temporal

        res = obj;

    elseif strcmp(obj.type,'until')

        % 1st rule in Lemma 4 in [1]
        if obj.from == 0 && obj.to == 0

            res = rewrite(obj.rhs);

        % 2nd rule in Table 1 in [1]
        elseif obj.from == 0

            I_ = interval(0,obj.to-1);
            phi1 = obj.lhs; phi2 = obj.rhs;

            res = rewrite(phi2) | (rewrite(phi1) & ...
                  rewrite(globally(phi1,I))) & (rewrite(finally(phi2,I) | ...
                  rewrite(next(until(phi1,phi2,I_),1))));

        % 1st rule in Table 1 in [1]
        else

            I_ = interval(obj.from-1,obj.to-1);
            phi1 = obj.lhs; phi2 = obj.rhs;

            res = rewrite(phi1) & rewrite(globally(phi1,I)) & ...
                            rewrite(next(until(phi1,phi2,I_),1));
        end

    elseif strcmp(obj.type,'release')

        % 2nd rule in Lemma 4 in [1]
        if obj.from == 0 && obj.to == 0

            res = rewrite(obj.rhs);

        % 4th rule in Table 1 in [1]
        elseif obj.from == 0

            I_ = interval(0,obj.to-1);
            phi1 = obj.lhs; phi2 = obj.rhs;

            res = rewrite(phi2) & (rewrite(phi1) | ...
                  rewrite(globally(phi2,I))) & (rewrite(finally(phi1,I) | ...
                  rewrite(next(release(phi1,phi2,I_),1))));

        % 3rd rule in Table 1 in [1]
        else

            I_ = interval(obj.from-1,obj.to-1);
            phi1 = obj.lhs; phi2 = obj.rhs;

            res = rewrite(phi1) | rewrite(finally(phi1,I)) | ...
                            rewrite(next(release(phi1,phi2,I_),1));

        end

    elseif strcmp(obj.type,'finally')

        % 3rd rule in Lemma 4 in [1]
        if strcmp(obj.lhs.type,'next')
    
            I_ = interval(obj.from,obj.to);

            res = next(rewrite(finally(obj.lhs.lhs,I_)),1);

        % 7th rule in Table 1 in [1]
        elseif strcmp(obj.lhs.type,'&') && obj.from == 0 && obj.to == 1

            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs;

            res = (rewrite(globally(phi1,I)) & rewrite(finally(phi2,I))) | ...
                  (rewrite(finally(phi1,I)) & rewrite(globally(phi2,I)));

        % 8th rule in Table 1 in [1]
        elseif strcmp(obj.lhs.type,'|') && obj.from == 0 && obj.to == 1

            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs;

            res = rewrite(finally(phi1,I)) | rewrite(finally(phi2,I));

        % 10th rule in Table 1 in [1]
        elseif strcmp(obj.lhs.type,'until') && obj.lhs.from == 0 && ...
                obj.lhs.to > 0 && obj.from == 0 && obj.to == 1

            I_ = interval(0,obj.lhs.to-1);
            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs;

            res = rewrite(finally(phi2,I)) | (rewrite(globally(phi1,I)) & ...
                    rewrite(next(phi2 | (phi1 & globally(phi1,I) & ...
                    finally(until(phi1,phi2,I_),I)),1)));

        % 9th rule in Table 1 in [1]
        elseif strcmp(obj.lhs.type,'until') && obj.lhs.from > 0 && ...
                obj.lhs.to > 0 && obj.from == 0 && obj.to == 1

           I_ = interval(obj.lhs.from-1,obj.lhs.to-1);
           phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs;

           res = rewrite(globally(phi1,I)) & rewrite(next(phi1 & ...
                    globally(phi1,I) & finally(until(phi1,phi2,I_),I),1));

        % 14th rule in Table 1 in [1]
        elseif strcmp(obj.lhs.type,'release') && obj.lhs.from == 0 && ...
                obj.lhs.to > 0 && obj.from == 0 && obj.to == 1

            I_ = interval(0,obj.lhs.to-1);
            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs;

            res = rewrite(finally(phi1 & phi2,I)) | (rewrite(globally(phi2,I)) & ...
                  rewrite(next(phi2 & (phi1 | (globally(phi2,I) & ...
                  finally(release(phi1,phi2,I_),I))),1)));

        % 13th rule in Table 1 in [1]
        elseif strcmp(obj.lhs.type,'release') && obj.lhs.from > 0 && ...
                obj.lhs.to > 0 && obj.from == 0 && obj.to == 1

            I_ = interval(obj.lhs.from-1,obj.lhs.to-1);
            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs;

            res = rewrite(finally(phi1,I)) | rewrite(next(phi1 | ...
                  finally(phi1,I) | finally(release(phi1,phi2,I_),I),1));  

        elseif obj.from ~= 0 || obj.to ~= 1

            res = rewrite(until(true,obj.lhs,interval(obj.from,obj.to)));

        elseif strcmp(obj.lhs.type,'finally') && obj.lhs.from == 0 && ...
                obj.lhs.to == 1 && obj.from == 0 && obj.to == 1

            I1 = interval(obj.from,obj.to);
            I2 = interval(obj.lhs.from,obj.lhs.to);

            res = rewrite(finally(until(true,obj.lhs.lhs,I2),I1));

        elseif strcmp(obj.lhs.type,'globally') && obj.lhs.from == 0 && ...
                obj.lhs.to == 1 && obj.from == 0 && obj.to == 1

            I1 = interval(obj.from,obj.to);
            I2 = interval(obj.lhs.from,obj.lhs.to);

            res = rewrite(finally(release(false,obj.lhs.lhs,I2),I1));

        elseif obj.lhs.temporal

            I_ = interval(obj.from,obj.to);

            res = rewrite(finally(rewrite(obj.lhs),I_));

        else
            res = obj;
        end

    elseif strcmp(obj.type,'globally')

        % 4th rule in Lemma 4 in [1]
        if strcmp(obj.lhs.type,'next')
    
            I_ = interval(obj.from,obj.to);

            res = next(rewrite(globally(obj.lhs.lhs,I_)),1);

        % 5th rule in Table 1 in [1]
        elseif strcmp(obj.lhs.type,'&') && obj.from == 0 && obj.to == 1

            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs;

            res = rewrite(globally(phi1,I)) & rewrite(globally(phi2,I));

        % 6th rule in Table 1 in [1]
        elseif strcmp(obj.lhs.type,'|') && obj.from == 0 && obj.to == 1

            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs;

            res = rewrite(globally(phi1,I)) | rewrite(globally(phi2,I));

        % 12th rule in Table 1 in [1]
        elseif strcmp(obj.lhs.type,'until') && obj.lhs.from == 0 && ...
                obj.lhs.to > 0 && obj.from == 0 && obj.to == 1

            I_ = interval(0,obj.lhs.to-1);
            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs;

            res = rewrite(globally(phi2,I)) | (rewrite(globally(phi1,I)) & ...
                  rewrite(next(phi2 | (phi1 & globally(phi1,I) & ...
                  globally(until(phi1,phi2,I_),I)),1)));

        % 11th rule in Table 1 in [1]
        elseif strcmp(obj.lhs.type,'until') && obj.lhs.from > 0 && ...
                obj.lhs.to > 0 && obj.from == 0 && obj.to == 1

            I_ = interval(obj.lhs.from-1,obj.lhs.to-1);
            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs;

            res = rewrite(globally(phi1,I)) & rewrite(next(phi1 & ...
                    globally(phi1,I) & globally(until(phi1,phi2,I_),I),1));

        % 16th rule in Table 1 in [1]
        elseif strcmp(obj.lhs.type,'release') && obj.lhs.from == 0 && ...
                obj.lhs.to > 0 && obj.from == 0 && obj.to == 1

            I_ = interval(0,obj.lhs.to-1);
            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs;

            res = rewrite(globally(phi2,I)) & (rewrite(globally(phi1,I) | ...
                  rewrite(next(phi2 & (phi1 | (globally(phi2,I) & ...
                  globally(release(phi1,phi2,I_),I))),1))));

        % 15th rule in Table 1 in [1]
        elseif strcmp(obj.lhs.type,'release') && obj.lhs.from > 0 && ...
                obj.lhs.to > 0 && obj.from == 0 && obj.to == 1

            I_ = interval(obj.lhs.from-1,obj.lhs.to-1);
            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs;

            res = rewrite(globally(phi1,I) | ...
                  rewrite(next(phi1 | globally(release(phi1,phi2,I_)))));
        
        elseif obj.from ~= 0 || obj.to ~= 1

            res = rewrite(release(false,obj.lhs,interval(obj.from,obj.to)));

        elseif strcmp(obj.lhs.type,'globally') && obj.lhs.from == 0 && ...
                obj.lhs.to == 1 && obj.from == 0 && obj.to == 1

            I1 = interval(obj.from,obj.to);
            I2 = interval(obj.lhs.from,obj.lhs.to);

            res = rewrite(globally(release(false,obj.lhs.lhs,I2),I1));

        elseif strcmp(obj.lhs.type,'finally') && obj.lhs.from == 0 && ...
                obj.lhs.to == 1 && obj.from == 0 && obj.to == 1

            I1 = interval(obj.from,obj.to);
            I2 = interval(obj.lhs.from,obj.lhs.to);

            res = rewrite(globally(until(true,obj.lhs.lhs,I2),I1));

        elseif obj.lhs.temporal

            I_ = interval(obj.from,obj.to);

            res = rewrite(globally(rewrite(obj.lhs),I_));

        else
            res = obj;
        end

    elseif strcmp(obj.type,'next')

        if strcmp(obj.lhs.type,'&')

            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs; time = obj.from;
            res = next(rewrite(phi1),time) & next(rewrite(phi2),time);

        elseif strcmp(obj.lhs.type,'|')

            phi1 = obj.lhs.lhs; phi2 = obj.lhs.rhs; time = obj.from;
            res = next(rewrite(phi1),time) | next(rewrite(phi2),time);

        else
            res = next(rewrite(obj.lhs),obj.from);
        end

    elseif strcmp(obj.type,'&')

        res = rewrite(obj.lhs) & rewrite(obj.rhs);

    elseif strcmp(obj.type,'|')

        res = rewrite(obj.lhs) | rewrite(obj.rhs);
    end
end

function res = scaleTime(obj,dt)
% recursive function to scale the time by the time step size dt

    scale = 1/dt;

    if ~obj.temporal

        res = obj;

    elseif ismember(obj.type,{'until','release'})

        res = obj;
        from = scale*res.from; to = scale*res.to;

        if abs(round(from) - from) > 1e-10 || abs(round(to) - to) > 1e-10
            throw(CORAerror( 'CORA:specialError', ...
                "Time step size does not divide times in STL formula!"));
        else
            res.from = int64(from); res.to = int64(to);
        end

        res.lhs = scaleTime(res.lhs,dt);
        res.rhs = scaleTime(res.rhs,dt);

    elseif ismember(obj.type,{'finally','globally'})

        res = obj;
        from = scale*res.from; to = scale*res.to;

        if abs(round(from) - from) > 1e-10 || abs(round(to) - to) > 1e-10
            throw(CORAerror( 'CORA:specialError', ...
                "Time step size does not divide times in STL formula!"));
        else
            res.from = int64(from); res.to = int64(to);
        end

        res.lhs = scaleTime(res.lhs,dt);

    elseif strcmp(obj.type,'next')

        res = obj;
        from = scale*res.from;

        if abs(round(from) - from) > 1e-10 
            throw(CORAerror( 'CORA:specialError', ...
                "Time step size does not divide times in STL formula!"));
        else
            res.from = int64(from);
        end

        res.lhs = scaleTime(res.lhs,dt);

    elseif strcmp(obj.type,'~')
        
        res = ~scaleTime(obj.lhs,dt);

    else

        res = obj;
        res.lhs = scaleTime(res.lhs,dt);
        res.rhs = scaleTime(res.rhs,dt);
    end
end

%------------- END OF CODE --------------