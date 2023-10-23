function [res,list] = stl2rtl(obj,dt,varargin)
% stl2rtl - convert STL formula to RTL formula according to Sec. 4.1 in [1] 
%
% Syntax:
%    [res,list] = stl2rtl(obj,dt)
%    [res,list] = stl2rtl(obj,dt,isSampled)
%
% Inputs:
%    obj - logic formula (class stl)
%    dt - time step size for sampling
%    isSampled - boolean specifying if formula is alread in sampled time
%
% Outputs:
%    res - resulting RTL formula (class rtl)
%    list - resulting RTL formula splitted into clauses
%
% Example: 
%    x = stl('x',2);
%    eq = ~(x(1) < 5 | globally(x(2) < 3,interval(0.1,0.2)));
%    eq_ = stl2rtl(eq,0.1)
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
% See also: stl, rtl

% Authors:       Niklas Kochdumper, Benedikt Seidl
% Written:       09-November-2022 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse input arguments
    isSampled = false;

    if nargin > 2 && ~isempty(varargin{1})
        isSampled = varargin{1};
    end

    % convert to sampled-time STL according to Sec. 4.2 in [1]
    if ~isSampled
        obj = sampledTime(obj,dt);
    end

    % bring to correct format required for Eq. (16) in [1]
    obj = combineNext(obj);
    obj = conjunctiveNormalForm(obj);

    % convert to RTL formula according to Lemma 2 in [1]
    clauses = getClauses(obj);
    res = []; list = {};

    for i = 1:length(clauses)
        
        [rtlTmp,listTmp] = aux_convertClause(clauses{i});

        res = res & rtlTmp;
        list = [list; {listTmp}];
    end
end


% Auxiliary functions -----------------------------------------------------

function [res,listRes] = aux_convertClause(obj)
% convert a single clause to RTL according to Lemma 2 in [1]

    % get all disjunctions
    clauses = getClauses(obj,'dnf');
    
    % group disjunctions with identical next operators together
    list = {};

    for i = 1:length(clauses)
        if strcmp(clauses{i}.type,'next')
            if length(list) < clauses{i}.from + 1
                list{clauses{i}.from+1} = {clauses{i}.lhs};
            else
                list{clauses{i}.from+1} = ...
                                [list{clauses{i}.from+1};{clauses{i}.lhs}];
            end
        else
            if length(list) < 1
                list{1} = clauses(i);
            else
                list{1} = [list{1};clauses(i)];
            end
        end
    end

    % loop over groups with identical next operator
    res = []; listRes = [];

    for i = 1:length(list)

        if ~isempty(list{i})

            % divide into finally, globally, and non-temporal
            listFin = {}; listGlob = {}; listNonTemp = {};
    
            for j = 1:length(list{i})
                if strcmp(list{i}{j}.type,'finally')
                    listFin{end+1} = list{i}{j}.lhs;
                elseif strcmp(list{i}{j}.type,'globally')
                    listGlob{end+1} = list{i}{j}.lhs;
                else
                    listNonTemp{end+1} = list{i}{j};
                end
            end

            % combine all finallies
            rho_ = [];

            for j = 1:length(listFin)
                rho_ = rho_ | listFin{j};
            end

            if ~isempty(rho_)
                tmp = next(rtl(rho_),i-0.5);
                res = res | tmp;
                listRes = [listRes; {tmp}];
            end

            % loop over all globallies
            for j = 1:length(listGlob)
                tmp = next(rtl(listGlob{j} | rho_),i-0.5);
                res = res | tmp;
                listRes = [listRes; {tmp}];
            end

            % loop over all non-temporal operators
            for j = 1:length(listNonTemp)
                tmp = next(rtl(listNonTemp{j}),i-1);
                res = res | tmp;
                listRes = [listRes; {tmp}];
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
