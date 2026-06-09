function [X0, spec, info] = priv_vnnlib2cora_v2(file)
% priv_vnnlib2cora_v2 - import specifications from VNN-LIB 2.0 .vnnlib files
%    Supports typed (declare-network ...) blocks, multiple inputs per network,
%    multi-dimensional bracket indexing X[i,j,k] (row-major), the == operator,
%    and multi-network specs with (isomorphic-to ...).
%
% Syntax:
%    [X0,spec,info] = priv_vnnlib2cora_v2(file)
%
% Inputs:
%    file - path to a .vnnlib file (VNN-LIB 2.0)
%
% Outputs:
%    X0   - cell array of input sets (interval for single-network box inputs,
%           polytope for multi-network or non-box inputs)
%    spec - specification on the (joint) output space
%    info - struct with .version, .networks, .inputDims, .outputDims,
%                       .totalIn, .totalOut, .isomorphism
%
% Reference:
%    - https://www.vnnlib.org/
%
% See also: vnnlib2cora, priv_vnnlib2cora_v1

% Authors:       Benedikt Kellner
% Written:       19-April-2026
% Last update:   08-May-2026 (BK, !=, declare-hidden, initialized, equal-to type, bounds checks)
%                22-May-2026 (BK, COO sparse accumulation — avoids dense C matrix for large inputs)
%                22-May-2026 (BK, regex fast path + vectorized interval extraction — speed optimization)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

text = fileread(file);
% strip line comments so the rest of the parser can ignore ';'
text = regexprep(text, ';[^\r\n]*', '');

% first pass: version header and (declare-network ...) blocks
version = aux_parseVersion(text);
[nets, varTable, totalIn, totalOut] = aux_parseNetworkDecls(text);
nNets = numel(nets);
if nNets == 0
    throw(CORAerror('CORA:converterIssue', ...
        sprintf('No networks declared in: %s', file)));
end

% accumulate constraint rows for each disjunct branch
data = aux_emptyData();

% fast path: bulk-extract all (assert (and (op VAR[idx] val) ...)) forms via regexp
[data, fastStarts, fastEnds] = aux_parseFastBoxConstraints(text, varTable, totalIn, totalOut, data);
for fi = 1:numel(fastStarts)
    text(fastStarts(fi):fastEnds(fi)) = ' ';
end

% second pass: walk remaining top-level forms (complex assertions only)
pos = 1;
while pos <= numel(text)
    pos = aux_skipWs(text, pos);
    if pos > numel(text) || text(pos) ~= '(', break; end

    % peek at the head keyword without consuming
    headPos = aux_skipWs(text, pos + 1);
    [head, ~] = aux_readToken(text, headPos);

    if strcmp(head, 'assert')
        % skip '( assert'
        pos = headPos + numel('assert');
        pos = aux_skipWs(text, pos);
        [pos, data] = aux_parseAssert(text, pos, data, varTable, totalIn, totalOut);
        pos = aux_skipWs(text, pos);
        pos = pos + 1; % closing ')'
    else
        % skip over any other top-level form by counting parentheses
        depth = 1;
        pos = pos + 1;
        while depth > 0 && pos <= numel(text)
            if text(pos) == '(', depth = depth + 1; end
            if text(pos) == ')', depth = depth - 1; end
            pos = pos + 1;
        end
    end
end

% build X0 from input polytopes
nIn = numel(data.polyInputs);
if nIn == 0
    throw(CORAerror('CORA:converterIssue', ...
        sprintf('No input constraints found in: %s', file)));
end
X0 = cell(1, nIn);
for b = 1:nIn
    br = data.polyInputs(b);
    C  = sparse(br.iRows, br.jCols, br.vals, br.nRows, totalIn);
    d  = br.dvec(:);
    if nNets == 1
        % single network: verify needs a box. An equality coupling two input
        % entries (X[i]==X[j], >=2 nonzero coeffs) is not axis-aligned -> reject.
        if br.nRowsEq > 0
            throw(CORAerror('CORA:converterIssue', 'Input set is not an interval.'));
        end
        X0{b} = aux_rowsToInterval(C, d, totalIn);
    else
        % multi-network: keep cross-variable couplings (X_f[i]==X_g[i]) as
        % native equality constraints (Ae*x==be) for the joint-network builder.
        Ce = sparse(br.iRowsEq, br.jColsEq, br.valsEq, br.nRowsEq, totalIn);
        X0{b} = polytope(full(C), d, full(Ce), br.dvecEq(:));
    end
end

% build spec from output polytopes (same reduction strategy as v1)
nOut = numel(data.polyOutputs);
if nOut == 0
    throw(CORAerror('CORA:converterIssue', ...
        sprintf('No output constraints found in: %s', file)));
end
Y = cell(1, nOut);
for b = 1:nOut
    br = data.polyOutputs(b);
    C  = full(sparse(br.iRows, br.jCols, br.vals, br.nRows, totalOut));
    d  = br.dvec(:);
    Y{b} = polytope(C, d);
end
if isscalar(Y)
    spec = specification(Y{1}, 'unsafeSet');
else
    Y_safe = safeSet2unsafeSet(Y);
    if length(Y_safe) < length(Y)
        Y = Y_safe; type = 'safeSet';
    else
        type = 'unsafeSet';
    end
    spec = [];
    for i = 1:length(Y)
        spec = add(spec, specification(Y{i}, type));
    end
end

% build info struct
inDims  = zeros(1, nNets);
outDims = zeros(1, nNets);
for k = 1:nNets
    inDims(k)  = sum(arrayfun(@(v) max(prod(v.shape), 1), nets(k).inputs));
    outDims(k) = max(prod(nets(k).output.shape), 1);
end
info.version     = version;
info.networks    = nets;
info.inputDims   = inDims;
info.outputDims  = outDims;
info.totalIn     = totalIn;
info.totalOut    = totalOut;
info.isomorphism = struct('target', {}, 'source', {}, 'type', {});
for k = 1:nNets
    if ~isempty(nets(k).equalTo)
        info.isomorphism(end+1) = struct( ...
            'target', nets(k).equalTo, 'source', nets(k).name, 'type', 'equal'); %#ok<AGROW>
    elseif ~isempty(nets(k).isomorphicTo)
        info.isomorphism(end+1) = struct( ...
            'target', nets(k).isomorphicTo, 'source', nets(k).name, 'type', 'isomorphic'); %#ok<AGROW>
    end
end

end


% Auxiliary functions -----------------------------------------------------

% --- First-pass helpers --------------------------------------------------

function version = aux_parseVersion(text)
% extract (vnnlib-version <X.Y>) header; default to '2.0' if absent
m = regexp(text, '\(\s*vnnlib-version\s*<\s*([\d\.]+)\s*>\s*\)', 'tokens', 'once');
if isempty(m)
    version = '2.0';
else
    version = m{1};
end
end

function [nets, varTable, totalIn, totalOut] = aux_parseNetworkDecls(text)
% parse every (declare-network ...) block and build variable lookup table

nets = repmat(aux_emptyNetwork(), 1, 0);
starts = strfind(text, '(declare-network');
for k = 1:numel(starts)
    % find matching closing ')' by counting paren depth
    depth = 0;
    for j = starts(k):numel(text)
        if text(j) == '(', depth = depth + 1; end
        if text(j) == ')', depth = depth - 1; end
        if depth == 0
            nets(end+1) = aux_parseOneNetwork(text(starts(k):j)); %#ok<AGROW>
            break;
        end
    end
end

netNames = {nets.name};
for k = 1:numel(nets)
    target = nets(k).isomorphicTo;
    if isempty(target), continue; end
    idx = find(strcmp(netNames, target), 1);
    if ~isempty(idx) && ~isempty(nets(idx).isomorphicTo)
        throw(CORAerror('CORA:converterIssue', ...
            sprintf('Chained network equivalence is not allowed: %s -> %s -> ...', ...
            nets(k).name, target)));
    end
end

% compute totalIn first so output offsets can reference it
totalIn = 0;
for k = 1:numel(nets)
    for j = 1:numel(nets(k).inputs)
        totalIn = totalIn + max(prod(nets(k).inputs(j).shape), 1);
    end
end

% register each input, output, and hidden variable with its joint-space column
varTable = containers.Map('KeyType','char','ValueType','any');
inOff = 0; outOff = 0;
for k = 1:numel(nets)
    for j = 1:numel(nets(k).inputs)
        inp = nets(k).inputs(j);
        if isKey(varTable, inp.name)
            throw(CORAerror('CORA:converterIssue', ...
                sprintf('Multiple declarations of variable: %s', inp.name)));
        end
        varTable(inp.name) = struct('kind','X','shape',inp.shape,'jointCol',inOff);
        inOff = inOff + max(prod(inp.shape), 1);
    end
    out = nets(k).output;
    if ~isempty(out.name)
        if isKey(varTable, out.name)
            throw(CORAerror('CORA:converterIssue', ...
                sprintf('Multiple declarations of variable: %s', out.name)));
        end
        varTable(out.name) = struct('kind','Y','shape',out.shape, ...
            'jointCol', totalIn + outOff);
        outOff = outOff + max(prod(out.shape), 1);
    end
    for j = 1:numel(nets(k).hidden)
        hid = nets(k).hidden(j);
        if isKey(varTable, hid.name)
            throw(CORAerror('CORA:converterIssue', ...
                sprintf('Multiple declarations of variable: %s', hid.name)));
        end
        varTable(hid.name) = struct('kind','H','shape',hid.shape,'jointCol',0);
    end
end
totalOut = outOff;
end

function net = aux_parseOneNetwork(block)
% extract name, inputs, output, and isomorphic-to from a block string

net = aux_emptyNetwork();

% network name: first identifier after 'declare-network'
m = regexp(block, '^\(declare-network\s+(\S+)', 'tokens', 'once');
if isempty(m)
    throw(CORAerror('CORA:converterIssue', 'Malformed declare-network block'));
end
net.name = m{1};

% all (declare-input NAME DTYPE [d1,d2,...] [initialized]) entries
inTokens = regexp(block, ...
    '\(declare-input\s+(\S+)\s+\S+\s+\[([^\]]*)\]\s*(initialized)?\s*\)', 'tokens');
emptyVar = struct('name', '', 'dtype', '', 'shape', [], 'initialized', false);
net.inputs = repmat(emptyVar, 1, 0);
for i = 1:numel(inTokens)
    v.name        = inTokens{i}{1};
    v.dtype       = 'real';
    v.shape       = aux_parseShape(strtrim(inTokens{i}{2}));
    v.initialized = ~isempty(inTokens{i}{3});
    net.inputs(end+1) = v; %#ok<AGROW>
end

% all (declare-hidden NAME DTYPE [shape] "onnxName") entries
hidTokens = regexp(block, ...
    '\(declare-hidden\s+(\S+)\s+\S+\s+\[([^\]]*)\]\s*"([^"]*)"\s*\)', 'tokens');
emptyHid = struct('name', '', 'dtype', 'real', 'shape', [], 'onnxName', '');
net.hidden = repmat(emptyHid, 1, 0);
for i = 1:numel(hidTokens)
    h.name     = hidTokens{i}{1};
    h.dtype    = 'real';
    h.shape    = aux_parseShape(strtrim(hidTokens{i}{2}));
    h.onnxName = hidTokens{i}{3};
    net.hidden(end+1) = h; %#ok<AGROW>
end

% (declare-output NAME DTYPE [SHAPE])
outTok = regexp(block, ...
    '\(declare-output\s+(\S+)\s+\S+\s+\[([^\]]*)\]\s*\)', 'tokens', 'once');
if ~isempty(outTok)
    net.output.name  = outTok{1};
    net.output.dtype = 'real';
    net.output.shape = aux_parseShape(strtrim(outTok{2}));
end

isoTok = regexp(block, '\(isomorphic-to\s+(\S+)\s*\)', 'tokens', 'once');
if ~isempty(isoTok)
    net.isomorphicTo = isoTok{1};
end

% equal-to also sets isomorphicTo for backward compatibility
eqTok = regexp(block, '\(equal-to\s+(\S+)\s*\)', 'tokens', 'once');
if ~isempty(eqTok)
    net.equalTo      = eqTok{1};
    net.isomorphicTo = eqTok{1};
end
end

function net = aux_emptyNetwork()
emptyVar    = struct('name', {}, 'dtype', {}, 'shape', {}, 'initialized', {});
emptyHidden = struct('name', {}, 'dtype', {}, 'shape', {}, 'onnxName', {});
net = struct('name', '', 'inputs', emptyVar, ...
    'output', struct('name','','dtype','','shape',[]), ...
    'isomorphicTo', '', ...
    'equalTo', '', ...
    'hidden', emptyHidden);
end

function data = aux_emptyData()
emptyBranch = aux_emptyBranch();
data = struct( ...
    'polyInputs',  emptyBranch, ...
    'polyOutputs', repmat(emptyBranch, 1, 0));  % output branches created lazily
end

function shape = aux_parseShape(str)
% parse comma-separated dimension list; empty string -> scalar shape 1
if isempty(str)
    shape = 1;
else
    shape = str2double(strsplit(str, ','));
end
end


% --- Second-pass helpers -------------------------------------------------

function [pos, data] = aux_parseAssert(text, pos, data, varTable, totalIn, totalOut)
% parse a parenthesized boolean expression at position pos

pos = aux_skip(text, pos, '(');
pos = aux_skipWs(text, pos);
posBeforeOp = pos;
[op, pos] = aux_readToken(text, pos);
pos = aux_skipWs(text, pos);

switch op
    case 'and'
        % all conjuncts share the current polytope branches
        while text(pos) ~= ')'
            [pos, data] = aux_parseAssert(text, pos, data, varTable, totalIn, totalOut);
            pos = aux_skipWs(text, pos);
        end

    case 'or'
        % each disjunct spawns independent polytope branches
        while text(pos) ~= ')'
            childData = aux_emptyData();
            [pos, childData] = aux_parseAssert(text, pos, childData, varTable, totalIn, totalOut);
            pos = aux_skipWs(text, pos);
            % append input branches only if the child added input constraints;
            % output-only disjuncts share the parent's input set
            if childData.polyInputs(1).nRows > 0 ...
                    || childData.polyInputs(1).nRowsEq > 0
                for j = 1:numel(childData.polyInputs)
                    data.polyInputs(end+1)  = childData.polyInputs(j);  %#ok<AGROW>
                end
            end
            for j = 1:numel(childData.polyOutputs)
                data.polyOutputs(end+1) = childData.polyOutputs(j); %#ok<AGROW>
            end
        end

    case '!='
        % != lhs rhs  ≡  (or (< lhs rhs) (> lhs rhs))
        % Split each existing branch into two rather than appending new empty ones.
        totalCols = totalIn + totalOut;
        [Clhs, dlhs, pos] = aux_parseExpr(text, pos, varTable, totalCols);
        pos = aux_skipWs(text, pos);
        [Crhs, drhs, pos] = aux_parseExpr(text, pos, varTable, totalCols);
        C = Clhs - Crhs;  d = drhs - dlhs;
        hasIn  = any(C(1:totalIn)     ~= 0);
        hasOut = any(C(totalIn+1:end) ~= 0);
        if hasIn && hasOut
            throw(CORAerror('CORA:notSupported', ...
                'Constraints mixing input and output variables are not supported'));
        end
        if hasIn || hasOut
            data = aux_splitBranchesOnNeq(data, C, d, hasIn, totalIn);
        end

    case {'<=', '<', '>=', '>', '=='}
        % linear comparison: parse both sides and add constraint rows
        totalCols = totalIn + totalOut;
        [Clhs, dlhs, pos] = aux_parseExpr(text, pos, varTable, totalCols);
        pos = aux_skipWs(text, pos);
        [Crhs, drhs, pos] = aux_parseExpr(text, pos, varTable, totalCols);

        % normalise to (Clhs - Crhs)*x <= drhs - dlhs
        % strict < and > are treated as non-strict (real-valued polytope representation)
        C = Clhs - Crhs;
        d = drhs  - dlhs;

        hasIn  = any(C(1:totalIn)       ~= 0);
        hasOut = any(C(totalIn+1:end)   ~= 0);
        if hasIn && hasOut
            throw(CORAerror('CORA:notSupported', ...
                'Constraints mixing input and output variables are not supported'));
        end

        if ~hasIn && ~hasOut
            % pure constant comparison — only check for contradiction
            if (any(strcmp(op, {'<=','<'})) && 0 > d) || ...
               (any(strcmp(op, {'>=','>'})) && 0 < -d) || ...
               (strcmp(op,'==') && d ~= 0)
                throw(CORAerror('CORA:converterIssue', ...
                    'Spec contains a contradictory constant constraint'));
            end
        else
            isInput = hasIn;
            switch op
                case {'<=', '<'}
                    data = aux_addRow(data,  C,  d, isInput, totalIn);
                case {'>=', '>'}
                    data = aux_addRow(data, -C, -d, isInput, totalIn);
                case '=='
                    if isInput && nnz(C(1:totalIn)) >= 2
                        % cross-variable input coupling (e.g. X_f[i]==X_g[j])
                        % -> keep as a native equality constraint (Ae*x==be).
                        data = aux_addEqRow(data, C, d, isInput, totalIn);
                    else
                        data = aux_addRow(data,  C,  d, isInput, totalIn);
                        data = aux_addRow(data, -C, -d, isInput, totalIn);
                    end
            end
        end

    otherwise
        % infix form: (lhs relop rhs) — backtrack and re-parse
        relOps = {'<=', '<', '>=', '>', '==', '!='};
        pos = posBeforeOp;
        totalCols = totalIn + totalOut;
        [Clhs, dlhs, pos] = aux_parseExpr(text, pos, varTable, totalCols);
        pos = aux_skipWs(text, pos);
        [op, pos] = aux_readToken(text, pos);
        pos = aux_skipWs(text, pos);
        if ~any(strcmp(op, relOps))
            throw(CORAerror('CORA:converterIssue', ...
                sprintf('Unknown assertion operator: %s', op)));
        end
        [Crhs, drhs, pos] = aux_parseExpr(text, pos, varTable, totalCols);
        % normalise and dispatch
        if strcmp(op, '!=')
            C = Clhs - Crhs;  d = drhs - dlhs;
            hasIn  = any(C(1:totalIn)     ~= 0);
            hasOut = any(C(totalIn+1:end) ~= 0);
            if hasIn && hasOut
                throw(CORAerror('CORA:notSupported', ...
                    'Constraints mixing input and output variables are not supported'));
            end
            if hasIn || hasOut
                data = aux_splitBranchesOnNeq(data, C, d, hasIn, totalIn);
            end
        else
            % relational operator: normalise to C*x <= d form and add rows
            C = Clhs - Crhs;  d = drhs - dlhs;
            hasIn  = any(C(1:totalIn)     ~= 0);
            hasOut = any(C(totalIn+1:end) ~= 0);
            if hasIn && hasOut
                throw(CORAerror('CORA:notSupported', ...
                    'Constraints mixing input and output variables are not supported'));
            end
            % constant or variable constraint
            if ~hasIn && ~hasOut
                if (any(strcmp(op, {'<=','<'})) && 0 > d) || ...
                   (any(strcmp(op, {'>=','>'})) && 0 < -d) || ...
                   (strcmp(op,'==') && d ~= 0)
                    throw(CORAerror('CORA:converterIssue', ...
                        'Spec contains a contradictory constant constraint'));
                end
            else
                % add halfspace row(s) to the appropriate polytope side
                isInput = hasIn;
                switch op
                    case {'<=', '<'}
                        data = aux_addRow(data,  C,  d, isInput, totalIn);
                    case {'>=', '>'}
                        data = aux_addRow(data, -C, -d, isInput, totalIn);
                    case '=='
                        if isInput && nnz(C(1:totalIn)) >= 2
                            % cross-variable input coupling -> native equality.
                            data = aux_addEqRow(data, C, d, isInput, totalIn);
                        else
                            data = aux_addRow(data,  C,  d, isInput, totalIn);
                            data = aux_addRow(data, -C, -d, isInput, totalIn);
                        end
                end
            end
        end
end

pos = aux_skipWs(text, pos);
pos = aux_skip(text, pos, ')');
end

function [C, d, pos] = aux_parseExpr(text, pos, varTable, totalCols)
% parse a linear arithmetic expression; returns row vector C and scalar d
% such that the expression equals C*[x;y] + d

pos = aux_skipWs(text, pos);
C = sparse(1, totalCols);
d = 0;

if text(pos) == '('
    % compound: (OP arg1 ...) — +, *, - are n-ary; unary (- x) is negation
    pos = pos + 1;
    pos = aux_skipWs(text, pos);
    [op, pos] = aux_readToken(text, pos);
    pos = aux_skipWs(text, pos);
    [C, d, pos] = aux_parseExpr(text, pos, varTable, totalCols);
    pos = aux_skipWs(text, pos);
    switch op
        case '+'
            while text(pos) ~= ')'
                [Ck, dk, pos] = aux_parseExpr(text, pos, varTable, totalCols);
                C = C + Ck; d = d + dk;
                pos = aux_skipWs(text, pos);
            end
        case '-'
            if text(pos) == ')'
                C = -C; d = -d; % unary negation
            else
                while text(pos) ~= ')'
                    [Ck, dk, pos] = aux_parseExpr(text, pos, varTable, totalCols);
                    C = C - Ck; d = d - dk;
                    pos = aux_skipWs(text, pos);
                end
            end
        case '*'
            while text(pos) ~= ')'
                [Ck, dk, pos] = aux_parseExpr(text, pos, varTable, totalCols);
                if ~any(C)
                    C = d * Ck; d = d * dk;
                elseif ~any(Ck)
                    C = C * dk; d = d * dk;
                else
                    throw(CORAerror('CORA:notSupported', ...
                        'Nonlinear (variable * variable) constraint is not supported'));
                end
                pos = aux_skipWs(text, pos);
            end
        otherwise
            throw(CORAerror('CORA:converterIssue', ...
                sprintf('Unknown arithmetic operator: %s', op)));
    end
    pos = aux_skip(text, pos, ')');
else
    % atom: number literal or variable reference
    [tok, pos] = aux_readToken(text, pos);
    c = tok(1);
    if (c >= '0' && c <= '9') || c == '-' || c == '.'
        val = str2double(tok);
        d = val;
    else
        % variable reference: NAME or NAME[i,j,k]
        indices = [];
        if pos <= numel(text) && text(pos) == '['
            pos = pos + 1; % skip '['
            [indices, pos] = aux_readIndexList(text, pos);
        end
        if ~isKey(varTable, tok)
            throw(CORAerror('CORA:converterIssue', ...
                sprintf('Unknown variable: %s', tok)));
        end
        vinfo = varTable(tok);
        if strcmp(vinfo.kind, 'H')
            throw(CORAerror('CORA:notSupported', ...
                sprintf('Constraints referencing hidden variable ''%s'' are not supported', tok)));
        end
        linIdx = aux_flatIdx(indices, vinfo.shape);
        C(vinfo.jointCol + linIdx) = 1;
    end
end
end

function data = aux_addRow(data, C, d, isInput, totalIn)
% append one constraint row to the COO accumulator for the appropriate side

if isInput
    field = 'polyInputs';
    row   = C(1:totalIn);
else
    field = 'polyOutputs';
    row   = C(totalIn+1:end);
end
if isempty(data.(field))
    data.(field)(1) = aux_emptyBranch();
end
% extract nonzeros from (sparse) row and append to every branch
nz  = find(row);
nNz = numel(nz);
v   = full(row(nz));
for i = 1:numel(data.(field))
    newR = data.(field)(i).nRows + 1;
    data.(field)(i).iRows(end+1:end+nNz) = newR;
    data.(field)(i).jCols(end+1:end+nNz) = nz;
    data.(field)(i).vals (end+1:end+nNz) = v;
    data.(field)(i).dvec (end+1)         = d;
    data.(field)(i).nRows                = newR;
end
end

function data = aux_addEqRow(data, C, d, isInput, totalIn)
% append one equality row (Ae*x==be) to every branch (cross-variable couplings)

if isInput
    field = 'polyInputs';
    row   = C(1:totalIn);
else
    field = 'polyOutputs';
    row   = C(totalIn+1:end);
end
if isempty(data.(field))
    data.(field)(1) = aux_emptyBranch();
end
nz  = find(row);
nNz = numel(nz);
v   = full(row(nz));
for i = 1:numel(data.(field))
    newR = data.(field)(i).nRowsEq + 1;
    data.(field)(i).iRowsEq(end+1:end+nNz) = newR;
    data.(field)(i).jColsEq(end+1:end+nNz) = nz;
    data.(field)(i).valsEq (end+1:end+nNz) = v;
    data.(field)(i).dvecEq (end+1)         = d;
    data.(field)(i).nRowsEq                = newR;
end
end

function b = aux_emptyBranch()
% inequality rows (.../nRows -> A*x<=b) and equality rows (...Eq -> Ae*x==be)
b = struct('iRows', [], 'jCols', [], 'vals', [], 'dvec', [], 'nRows', 0, ...
    'iRowsEq', [], 'jColsEq', [], 'valsEq', [], 'dvecEq', [], 'nRowsEq', 0);
end


% --- Character-level utilities -------------------------------------------

function pos = aux_skipWs(text, pos)
% advance past whitespace characters
while pos <= numel(text) && ...
        (text(pos) == ' ' || text(pos) == char(9) || ...
         text(pos) == char(10) || text(pos) == char(13))
    pos = pos + 1;
end
end

function pos = aux_skip(text, pos, ch)
% assert the character at pos equals ch and advance
if pos > numel(text) || text(pos) ~= ch
    throw(CORAerror('CORA:converterIssue', ...
        sprintf('Expected ''%s'' at position %d (got ''%s'')', ...
        ch, pos, text(min(pos, numel(text))))));
end
pos = pos + 1;
end

function [tok, pos] = aux_readToken(text, pos)
% read one token (non-whitespace, non-delimiter run); < and > excluded so
% that <=, >= are read as single tokens
delim = '()[],; ';
start = pos;
while pos <= numel(text)
    c = text(pos);
    if any(c == delim) || c == char(9) || c == char(10) || c == char(13)
        break;
    end
    pos = pos + 1;
end
tok = text(start:pos-1);
end

function [indices, pos] = aux_readIndexList(text, pos)
% read comma-separated integers up to and including the closing ']'
indices = [];
while pos <= numel(text) && text(pos) ~= ']'
    if text(pos) == ','
        pos = pos + 1;
        continue;
    end
    [tok, pos] = aux_readToken(text, pos);
    indices(end+1) = str2double(tok); %#ok<AGROW>
end
pos = pos + 1; % skip ']'
end

function linIdx = aux_flatIdx(indices, shape)
% row-major (C order) flat index; VNN-LIB is 0-based, MATLAB 1-based
if isempty(indices)
    linIdx = 1;
    return;
end
if numel(indices) ~= numel(shape)
    % scalar variable indexed as [0]
    if isscalar(shape) && isscalar(indices) && indices == 0
        linIdx = 1;
        return;
    end
    throw(CORAerror('CORA:converterIssue', ...
        sprintf('Index dim %d does not match shape dim %d', ...
        numel(indices), numel(shape))));
end
linIdx = 0;
for k = 1:numel(indices)
    if indices(k) < 0 || indices(k) >= shape(k)
        throw(CORAerror('CORA:converterIssue', ...
            sprintf('Index out of bounds: index %d (0-based) in dimension %d with size %d', ...
            indices(k), k, shape(k))));
    end
    linIdx = linIdx + indices(k) * prod(shape(k+1:end));
end
linIdx = linIdx + 1;
end

function I = aux_rowsToInterval(C, d, n)
% Fast path: if all rows are axis-aligned ±1 half-spaces, extract bounds
% directly without building a polytope.  Falls back to polytope+representsa_
% for files that include redundant non-axis-aligned constraints.
%
% Vectorized: find all nonzeros at once to avoid row-by-row O(n^2) scan.
nzPerRow = sum(C ~= 0, 2);
if any(nzPerRow ~= 1) || any(abs(C(C ~= 0)) ~= 1)
    % non-axis-aligned constraint present — fall back
    [isBox, I] = representsa_(polytope(C, d), 'interval', eps);
    if ~isBox
        throw(CORAerror('CORA:converterIssue', 'Input set is not an interval.'));
    end
    return;
end

% All rows are ±1 axis-aligned: extract lb/ub directly.
lb = -inf(n, 1);
ub =  inf(n, 1);
[rows, cols] = find(C);  % one vectorized scan
pos = C(sub2ind(size(C), rows, cols)) > 0;  % true -> upper bound row
ubMask = logical(pos);
lbMask = ~ubMask;
if any(ubMask)
    ub = min(ub, accumarray(cols(ubMask), d(rows(ubMask)),   [n,1], @min,  inf));
end
if any(lbMask)
    lb = max(lb, accumarray(cols(lbMask), -d(rows(lbMask)), [n,1], @max, -inf));
end
I = interval(lb, ub);
end

function data = aux_splitBranchesOnNeq(data, C, d, isInput, totalIn)
% Splits each branch in two (C*x<=d and -C*x<=-d) rather than appending new
% empty ones — avoids the stray initial branch that the 'or' pattern produces.
if isInput
    field = 'polyInputs';
    row = C(1:totalIn);
else
    field = 'polyOutputs';
    row = C(totalIn+1:end);
end
nz  = find(row);
nNz = numel(nz);
v   = full(row(nz));
% clone current branches and append +row / -row to each side
branchesA = data.(field);
branchesB = data.(field);
for b = 1:numel(branchesA)
    newR = branchesA(b).nRows + 1;
    branchesA(b).iRows(end+1:end+nNz) = newR;
    branchesA(b).jCols(end+1:end+nNz) = nz;
    branchesA(b).vals (end+1:end+nNz) =  v;
    branchesA(b).dvec (end+1)         =  d;
    branchesA(b).nRows                = newR;
    branchesB(b).iRows(end+1:end+nNz) = newR;
    branchesB(b).jCols(end+1:end+nNz) = nz;
    branchesB(b).vals (end+1:end+nNz) = -v;
    branchesB(b).dvec (end+1)         = -d;
    branchesB(b).nRows                = newR;
end
data.(field) = [branchesA, branchesB];
end

function [data, matchStarts, matchEnds] = ...
        aux_parseFastBoxConstraints(text, varTable, totalIn, totalOut, data)
% Bulk-extract (assert (and (op NAME[idx] val) (op NAME[idx] val))) patterns
% via regexp and append them directly to the COO accumulators.
% Fully vectorized: no MATLAB interpreted loop over assertions.
% Returns matched position ranges so the caller can blank them in text.

matchStarts = [];
matchEnds   = [];

numPat = '[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?';
opPat  = '[<>]=?|==';
idxPat = '[0-9, ]+';
onePat = ['(' opPat ')\s+(\w+)\[(' idxPat ')\]\s+(' numPat ')'];
fullPat = ['\(assert\s+\(and\s+\(' onePat '\)\s+\(' onePat '\)\s*\)\s*\)'];

[toks, matchStarts, matchEnds] = regexp(text, fullPat, 'tokens', 'start', 'end');
if isempty(toks)
    return;
end

% Flatten: toks{k} = {op1,name1,idx1,val1, op2,name2,idx2,val2}
% Stack into 8-column cell array, then split sides and interleave.
all8 = vertcat(toks{:});          % N × 8 cell
% Each row: [op1 name1 idx1 val1  op2 name2 idx2 val2]
% Interleave both sides so we process them together.
allOps   = [all8(:,1); all8(:,5)];   % 2N × 1
allNames = [all8(:,2); all8(:,6)];
allIdxs  = [all8(:,3); all8(:,7)];
allVals  = str2double([all8(:,4); all8(:,8)]);  % 2N × 1 (vectorized)

% Process per unique variable (typically just 2-3 variables).
uniqueVars = unique(allNames);
for v = 1:numel(uniqueVars)
    vname = uniqueVars{v};
    if ~isKey(varTable, vname), continue; end
    vinfo = varTable(vname);
    if strcmp(vinfo.kind, 'H'), continue; end

    mask  = strcmp(allNames, vname);
    ops   = allOps(mask);
    vals  = allVals(mask);
    idxStrs = allIdxs(mask);
    M     = sum(mask);

    % Vectorized index computation via sscanf over joined string.
    ndim    = numel(vinfo.shape);
    joined  = strjoin(idxStrs, ',');
    rawNums = sscanf(joined, '%f,');           % ndim*M values
    idxMat  = reshape(rawNums, ndim, M)';      % M × ndim (0-based)
    strides = [fliplr(cumprod(fliplr(vinfo.shape(2:end)))), 1];
    linIdxs = idxMat * strides(:) + 1;        % M × 1, 1-based
    cols    = vinfo.jointCol + linIdxs;        % M × 1

    % Operator → sign/rhs transformation (vectorized).
    isLE = strcmp(ops, '<=') | strcmp(ops, '<');
    isGE = strcmp(ops, '>=') | strcmp(ops, '>');
    isEQ = ~isLE & ~isGE;

    sgns = zeros(M, 1);
    sgns(isLE) =  1;
    sgns(isGE) = -1;
    dVals = vals;
    dVals(isGE) = -vals(isGE);

    % Split into input vs. output.
    isInp = cols <= totalIn;

    % EQ constraints produce two rows each; handle separately (rare).
    if any(isEQ)
        eqCols  = cols(isEQ);
        eqDVals = vals(isEQ);
        eqInp   = isInp(isEQ);
        data = aux_appendEqRows(data, eqCols, eqDVals, eqInp, totalIn);
        % Exclude EQ rows from the main append below.
        isLE(isEQ) = false;
        isGE(isEQ) = false;
        isInp(isEQ) = false;
    end

    % Append all <= and >= rows for this variable in one shot.
    for fldIdx = 1:2
        if fldIdx == 1
            field   = 'polyInputs';
            sel     = isInp;
            fldCols = cols(sel);
        else
            field   = 'polyOutputs';
            sel     = ~isInp & (isLE | isGE);
            fldCols = cols(sel) - totalIn;
        end
        if ~any(sel), continue; end
        if isempty(data.(field))
            data.(field)(1) = aux_emptyBranch();
        end
        % bulk-append all rows for this variable/side in one assignment
        fldSgns  = sgns(sel);
        fldDVals = dVals(sel);
        nNew = sum(sel);
        for i = 1:numel(data.(field))
            startR = data.(field)(i).nRows;
            data.(field)(i).iRows(end+1:end+nNew) = (startR+1):(startR+nNew);
            data.(field)(i).jCols(end+1:end+nNew) = fldCols(:)';
            data.(field)(i).vals (end+1:end+nNew) = fldSgns(:)';
            data.(field)(i).dvec (end+1:end+nNew) = fldDVals(:)';
            data.(field)(i).nRows = startR + nNew;
        end
    end
end
end

function data = aux_appendEqRows(data, cols, dVals, isInp, totalIn)
% Append == constraints (two halfspaces each) — called only for rare == ops.
for k = 1:numel(cols)
    col = cols(k);
    val = dVals(k);
    if isInp(k)
        field    = 'polyInputs';
        localCol = col;
    else
        field    = 'polyOutputs';
        localCol = col - totalIn;
    end
    if isempty(data.(field))
        data.(field)(1) = aux_emptyBranch();
    end
    for i = 1:numel(data.(field))
        newR = data.(field)(i).nRows + 1;
        data.(field)(i).iRows(end+1:end+2) = [newR newR+1];
        data.(field)(i).jCols(end+1:end+2) = [localCol localCol];
        data.(field)(i).vals (end+1:end+2) = [1 -1];
        data.(field)(i).dvec (end+1:end+2) = [val -val];
        data.(field)(i).nRows = newR + 1;
    end
end
end

% ------------------------------ END OF CODE ------------------------------
