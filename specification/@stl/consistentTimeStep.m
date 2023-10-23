function res = consistentTimeStep(obj,dt)
% consistentTimeStep - make STL formula consistent with the given time step 
%
% Syntax:
%    res = consistentTimeStep(obj,dt)
%
% Inputs:
%    obj - logic formula (class stl)
%    dt - time step size
%
% Outputs:
%    res - resulting consistent logic formula (class stl)
%
% Example: 
%    x = stl('x',2);
%    eq = globally(x(1) + x(2) < 3, interval(1.5, 3.2));
%
%    res = consistentTimeStep(eq,1)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl, stl/sampledTime

% Authors:       Niklas Kochdumper
% Written:       09-November-2022 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    scale = 1/dt;

    if ~obj.temporal

        res = obj;

    elseif strcmp(obj.type,'until')

        res = obj;
        from = scale*res.from; to = scale*res.to;

        if abs(round(from) - from) > 1e-10 || abs(round(to) - to) > 1e-10
            res.from = ceil(from)*dt; res.to = floor(to)*dt;
        else
            res.from = round(from)*dt; res.to = round(to)*dt;
        end

        res.lhs = consistentTimeStep(res.lhs,dt);
        res.rhs = consistentTimeStep(res.rhs,dt);

    elseif strcmp(obj.type,'release')

        res = obj;
        from = scale*res.from; to = scale*res.to;

        if abs(round(from) - from) > 1e-10 || abs(round(to) - to) > 1e-10
            res.from = floor(from)*dt; res.to = ceil(to)*dt;
        else
            res.from = round(from)*dt; res.to = round(to)*dt;
        end

        res.lhs = consistentTimeStep(res.lhs,dt);
        res.rhs = consistentTimeStep(res.rhs,dt);

    elseif strcmp(obj.type,'globally')

        res = obj;
        from = scale*res.from; to = scale*res.to;

        if abs(round(from) - from) > 1e-10 || abs(round(to) - to) > 1e-10
            res.from = floor(from)*dt; res.to = ceil(to)*dt;
        else
            res.from = round(from)*dt; res.to = round(to)*dt;
        end

        res.lhs = consistentTimeStep(res.lhs,dt);

    elseif strcmp(obj.type,'finally')

        res = obj;
        from = scale*res.from; to = scale*res.to;

        if abs(round(from) - from) > 1e-10 || abs(round(to) - to) > 1e-10
            res.from = ceil(from)*dt; res.to = floor(to)*dt;
        else
            res.from = round(from)*dt; res.to = round(to)*dt;
        end

        res.lhs = consistentTimeStep(res.lhs,dt);

    elseif strcmp(obj.type,'next')

        res = obj;
        from = scale*res.from;

        if abs(round(from) - from) > 1e-10 
            inner = consistentTimeStep(res.lhs,dt);
            res = globally(inner,interval(floor(from)*dt,ceil(from)*dt));
        else
            res.from = round(from)*dt;
            res.lhs = consistentTimeStep(res.lhs,dt);
        end

    elseif strcmp(obj.type,'~')
        
        res = ~consistentTimeStep(obj.lhs,dt);

    else

        res = obj;
        res.lhs = consistentTimeStep(res.lhs,dt);
        res.rhs = consistentTimeStep(res.rhs,dt);
    end
end

% ------------------------------ END OF CODE ------------------------------
