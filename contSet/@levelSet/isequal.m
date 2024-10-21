function res = isequal(ls,S,varargin)
% isequal - checks if a levelSet objects is equal to another set or point
%
% Syntax:
%    res = isequal(ls,S)
%    res = isequal(ls,S,tol)
%
% Inputs:
%    ls - levelSet object
%    S - contSet object, numeric
%    tol - tolerance (optional)
%
% Outputs:
%    res - true/false
%
% Example:
%    % init symbolic variables
%    syms a b c x y z
% 
%    % init level set
%    eq = -x^2 - y^2 + 5;
%    ls1 = levelSet(eq,[x;y],'<=');
% 
%    % different variable names
%    eq = -a^2 - b^2 + 5;
%    ls1_ = levelSet(eq,[a;b],'<=');
% 
%    % compare equal sets
%    isequal(ls1,ls1_)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       11-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(2,3);

% default values
tol = setDefaultValues({eps},varargin);

% check input arguments
inputArgsCheck({{ls,'att',{'levelSet','numeric'}};
                {S,'att',{'contSet','numeric'}};
                {tol,'att','numeric',{'scalar','nonnegative','nonnan'}}});

% ensure that numeric is second input argument
[ls,S] = reorderNumeric(ls,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < ls.precedence
    res = isequal(S,ls,tol);
    return
end

% ambient dimensions must match
if ~equalDimCheck(ls,S,true)
    res = false;
    return
end

% levelSet-levelSet case
if isa(S,'levelSet')
    res = aux_isequal_levelSet(ls,S,tol);
    return
end

throw(CORAerror('CORA:noops',ls,S));

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isequal_levelSet(ls,S,tol)

% reduce both to minimal representation
ls = compact_(ls,'all',eps);
S = compact_(S,'all',eps);

% check if number of equations are equal
if length(ls.eq) ~= length(S.eq)
    res = false; return
end

% check if number of used (!) variables is equal
usedVarsTotal = length(aux_usedVariables(ls));
if usedVarsTotal ~= length(aux_usedVariables(ls))
    res = false; return
end

% check if symbolic expressions and comparison operator match
idxInLS2 = false(length(ls.eq));

% init map for symbolic variable name comparison:
% - (:,1): first level set
% - (:,2): second level set
% - (:,3): chosen variable names (enabling comparison)
symNames = [sym('notassigned_ls1_',[usedVarsTotal,1],'real'), ...
    sym('notassigned_ls2_',[usedVarsTotal,1],'real'), ...
    sym('xyz_',[usedVarsTotal,1],'real')];
% logical arrays which of the substituted variable names have been assigned
assigned = false(usedVarsTotal,1);


% loop over all equations from first level set
for i=1:length(ls.eq)

    % assume: no match found
    found = false;

    % number of variables used in i-th equation of first level set
    usedVars1 = symvar(ls.eq(i)).';
    nrUsedVars1 = length(usedVars1);

    % find out which variable names should be substituted for given
    % symbolic variables (keep correspondance over mulitple equations!)
    if i == 1
        % fill up from start (without loss of generality)
        symNames(1:nrUsedVars1,1) = usedVars1;
        % names of substituted variables
        subNames = symNames(1:length(usedVars1),3);
        % save which ones have been assigned
        assigned(1:nrUsedVars1,1) = true;
        % indices
        idx = 1:nrUsedVars1;

    else
        % check if any symbolic variables have already been used

        % names of substituted variables
        subNames = [];
        % indices of substituted variables
        idx = zeros(nrUsedVars1,1);

        for v=1:nrUsedVars1
            idx_temp = logical(symNames(:,1) - usedVars1(v) == 0);

            if any(idx_temp)
                % variable has already been assigned before -> read out
                subNames = [subNames; symNames(idx_temp,3)];

                % index for substitution
                idx(v) = find(idx_temp,1,'first');
            else
                % variable has not been assigned -> pick next available
                % variable (without loss of generality)
                idx_free = find(~assigned,1,'first');
                assigned(idx_free) = true;
                subNames = [subNames; symNames(idx_free,3)];

                % assign variable
                symNames(idx_free,1) = usedVars1(v);

                % index for substitution
                idx(v) = idx_free;
            end
        end

    end

    % substitute variable names for comparison in first level set
    ls1_subbed = subs(ls.eq(i),usedVars1,subNames);
    

    % loop over all equations from second level set
    for j=1:length(S.eq)
        % skip if j-th entry has already been matched to an equation of the
        % first level set

        if ~idxInLS2(j)
            % j-th equation of second level set not yet matched

            % check if same comparison operator used (different handling
            % depending on whether compOp is a cell-array)
            if ~iscell(ls.compOp) && ~strcmp(ls.compOp,S.compOp) ...
                    || iscell(ls.compOp) && ~strcmp(ls.compOp{i},S.compOp{j})
                % try next equation of second level set
                continue;
            end

            % variable names used in j-th equation
            usedVars2 = symvar(S.eq(j)).';
            % number of symbolic variables in j-th equation
            nrUsedVars2 = length(usedVars2);

            % check if same number of symbolic variables used
            if nrUsedVars2 ~= nrUsedVars1
                % try next equation of second level set
                continue;
            end

            % substitute same variable names in second level set's equation
            ls2_subbed = subs(S.eq(j),usedVars2,subNames);

            % check if difference between two equations equals 0
            if logical(simplify(ls1_subbed - ls2_subbed) == 0)
                % match has been found
                found = true; idxInLS2(j) = true;
                % save matched variable names for subsequent comparisons
                symNames(idx,2) = usedVars2;
                break
            end
        end
    end

    if ~found
        % i-th equation in first level set has no match in second level set
        res = false; return
    end
end

% all checks ok
res = true;

end

function usedVars = aux_usedVariables(ls)
% returns number of variables that are actually used in the level set
% equations

usedVars = [];
% loop over all equations
for i=1:length(ls.eq)
    usedVars = [usedVars,symvar(ls.eq(i))];
end
% remove duplicates
usedVars = unique(usedVars);

end

% ------------------------------ END OF CODE ------------------------------
