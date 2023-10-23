function res = in(obj,S)
% in - set containment operator for Signal Temporal Logic
%
% Syntax:
%    res = in(obj,S)
%
% Inputs:
%    obj - logic formula (class stl)
%    S - set (class interval, zonotope, polytope, conZonotope,
%               zonoBundle, halfspace, or polygon)
%
% Outputs:
%    res - resulting stl formula (class stl)
%
% Example: 
%    x = stl('x',2);
%    pgon = polygon.generateRandom();
%    eq = finally(in(x,pgon),interval(0.1,0.2))
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

% check input arguments
if ~isa(obj,'stl') || ~ismember(obj.type,{'variable','concat'})
    throw(CORAerror('CORA:notSupported',...
                  'This operation is not supported for stl objects!'));
end

if (strcmp(obj.type,'variable') && dim(S) ~= length(obj.variables)) || ...
        strcmp(obj.type,'concat') && dim(S) ~= length(obj.var)
    throw(CORAerror('CORA:wrongValue','second', ...
                  'dimensions of set and stl object have to match"!'));
end

% conversions for different types of sets
if isa(S,'polygon')

    list = splitIntoConvexSets(S);

    res = [];
    for i = 1:length(list)
        res = res | polytope2stl(obj,polytope(list{i}));
    end

elseif isa(S,'interval') || isa(S,'zonotope') || ...
       isa(S,'polytope') || isa(S,'conZonotope') || ...
       isa(S,'zonoBundle') || isa(S,'halfspace')

    res = polytope2stl(obj,polytope(S));

elseif isa(S,'levelSet') || isa(S,'ellipsoid')

    S = levelSet(S);

    if size(S.eq,1) == 1
        res = aux_inequalityLevelSet(S.funHan(obj),S.compOp);
    else
        for i = 1:size(S.eq,1)
            eq = matlabFunction(S.eq(i,1),'Vars',{S.vars});
            if i == 1
                res = aux_inequalityLevelSet(eq(obj),S.compOp{i});
            else
                res = res & aux_inequalityLevelSet(eq(obj),S.compOp{i});
            end
        end
    end

elseif isa(S,'capsule')

    n = dim(S);
    x = sym('x',[n,1]);

    % cylinder of the capsule
    x_ = (x - S.c) - S.g * (S.g'*(x - S.c)) / norm(S.g)^2;
    eq1 = sum(x_.^2) - S.r^2;

    eq2 = S.g' * x - S.g'*(S.c + S.g);
    eq3 = -S.g' * x + S.g'*(S.c - S.g);

    ls1 = levelSet([eq1;eq2;eq3],x,{'<=';'<=';'<='});

    % the two spheres of the capsule
    ls2 = levelSet(ellipsoid(eye(n)*S.r^2) + S.c + S.g);
    ls3 = levelSet(ellipsoid(eye(n)*S.r^2) + S.c - S.g);

    % construct logic equations
    res = in(obj,ls1) | in(obj,ls2) | in(obj,ls3);

else
    throw(CORAerror('CORA:notSupported',...
         'This operation is not supported for this type of set representation!'));
end


% Auxiliary functions -----------------------------------------------------

function res = aux_inequalityLevelSet(obj,op)
% construct an STL formula for the inequality of a level set constraint

    if strcmp(op,'<')
        res = obj < 0;
    elseif strcmp(op,'<=')
        res = obj <= 0;
    elseif strcmp(op,'>')
        res = obj > 0;
    elseif strcmp(op,'>=')
        res = obj >= 0;
    else
        throw(CORAerror('CORA:notSupported',...
         'Level sets with equality constraints cannot be converted to STL!'));
    end

% ------------------------------ END OF CODE ------------------------------
