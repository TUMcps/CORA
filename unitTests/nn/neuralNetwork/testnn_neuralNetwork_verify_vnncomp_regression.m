function res = testnn_neuralNetwork_verify_vnncomp_regression()
% testnn_neuralNetwork_verify_vnncomp_regression - test writeCounterexample
%    for VNN-LIB 1.x and 2.0 output (format, round-trip, multi-network layout,
%    v1 fallback, save/load), using vnnlib2cora fixtures
%
% Syntax:
%    res = testnn_neuralNetwork_verify_vnncomp_regression()
%
% Inputs:
%    -
%
% Outputs:
%    res - true if all tests pass, false otherwise
%
% Other m-files required: writeCounterexample, vnnlib2cora
% Subfunctions: aux_write, aux_parseV1, aux_parseV2, aux_headerOrder
% MAT-files required: none
%
% See also: writeCounterexample, run_instance,
%    testnn_neuralNetwork_verify_vnncomp_ce_format

% Authors:       Benedikt Kellner
% Written:       05-June-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% temp workspace
testDir = tempname;
mkdir(testDir);
cleanupDir = onCleanup(@() rmdir(testDir, 's'));

% writeCounterexample lives in the vnncomp directory
vnncompPath = fullfile(CORAROOT, 'examples', 'nn', 'vnncomp');
addpath(vnncompPath);

v2Dir = fullfile(CORAROOT, 'models', 'Cora', 'nn', 'unitTests', 'vnnlib', 'v2');

% ----- v1: empty vnnlibInfo falls back to legacy format ------------------
x_ = [0.1; 0.2; 0.3]; y_ = [1.0; 2.0];
content = aux_write(testDir, x_, y_, []);
[xp, yp] = aux_parseV1(content);
assert(startsWith(content, ['sat' newline]), 'v1: missing sat header');
assert(contains(content, sprintf('%s(', newline)), 'v1: missing outer paren');
assert(isequal(xp, x_) && isequal(yp, y_), 'v1: empty-info round-trip failed');

% ----- v1: struct without a version field also falls back ----------------
content = aux_write(testDir, x_, y_, struct());
assert(~isempty(regexp(content, '\(X_0\s', 'once')), ...
    'v1: missing-version-field did not fall back to legacy format');

% ----- v1: version ~= 2.0 uses legacy format -----------------------------
content = aux_write(testDir, x_, y_, struct('version', '1.0'));
assert(~isempty(regexp(content, '\(X_0\s', 'once')), ...
    'v1: version 1.0 did not use legacy format');

% ----- v1: large vectors keep every element (%f is 6-decimal lossy) ------
xL = rand(100, 1); yL = rand(50, 1);
[xp, yp] = aux_parseV1(aux_write(testDir, xL, yL, []));
assert(numel(xp) == 100 && numel(yp) == 50, 'v1: large element count wrong');
assert(max(abs(xp - xL)) <= 1e-6 && max(abs(yp - yL)) <= 1e-6, ...
    'v1: large round-trip outside %f precision');

% ----- v2: single-network metadata from vnnlib2cora ----------------------
% multidim_basic declares X real [2,3] and Y real [2].
f = fullfile(v2Dir, 'multidim_basic.vnnlib');
assert(isfile(f), 'fixture multidim_basic.vnnlib missing');
[~, ~, info] = vnnlib2cora(f);
x_ = (1:6)' / 10; y_ = [3.5; -2.0];
content = aux_write(testDir, x_, y_, info);
% header lines must show the declared dtype ('real') and the shape
assert(contains(content, 'X real [2,3]'), 'v2: wrong input header');
assert(contains(content, 'Y real [2]'), 'v2: wrong output header');
assert(isempty(regexp(content, '\(X_\d', 'once')), 'v2: leaked v1 syntax');
[xp, yp] = aux_parseV2(content, {'X'}, {'Y'});
assert(isequal(xp, x_) && isequal(yp, y_), 'v2: single-network round-trip failed');

% ----- v2: multi-network metadata (joint layout + decl order) ------------
% equal_to declares f(Xf[3]->Yf[2]) then g(Xg[3]->Yg[2]).
f = fullfile(v2Dir, 'equal_to.vnnlib');
assert(isfile(f), 'fixture equal_to.vnnlib missing');
[~, ~, info] = vnnlib2cora(f);
assert(numel(info.networks) == 2, 'equal_to should declare 2 networks');
% joint layout: x_ = [Xf(3); Xg(3)], y_ = [Yf(2); Yg(2)]
xf = [0.11; 0.12; 0.13]; xg = [0.21; 0.22; 0.23];
yf = [1.1; 1.2];          yg = [2.1; 2.2];
x_ = [xf; xg]; y_ = [yf; yg];
content = aux_write(testDir, x_, y_, info);

% all four declared variables present with their names/shapes
nf = info.networks(1); ng = info.networks(2);
assert(contains(content, sprintf('%s real [3]', nf.inputs(1).name)), 'v2 multinet: Xf header');
assert(contains(content, sprintf('%s real [2]', nf.output.name)),    'v2 multinet: Yf header');
assert(contains(content, sprintf('%s real [3]', ng.inputs(1).name)), 'v2 multinet: Xg header');
assert(contains(content, sprintf('%s real [2]', ng.output.name)),    'v2 multinet: Yg header');

% declaration order must be Xf, Yf, Xg, Yg (network f fully, then g)
order = aux_headerOrder(content);
assert(isequal(order, {nf.inputs(1).name, nf.output.name, ...
    ng.inputs(1).name, ng.output.name}), ...
    'v2 multinet: header order does not follow declaration order');

% values round-trip per variable against the joint layout
[xp, yp] = aux_parseV2(content, ...
    {nf.inputs(1).name, ng.inputs(1).name}, {nf.output.name, ng.output.name});
assert(isequal(xp, x_), 'v2 multinet: input round-trip failed');
assert(isequal(yp, y_), 'v2 multinet: output round-trip failed');

% ----- backward compat: old .mat without vnnlibInfo loads as [] ----------
matFile = fullfile(testDir, 'old.mat');
foo = 1; save(matFile, 'foo'); %#ok<NASGU>
clear vnnlibInfo;
warning('off', 'MATLAB:load:variableNotFound');
load(matFile, 'vnnlibInfo');
warning('on', 'MATLAB:load:variableNotFound');
if ~exist('vnnlibInfo', 'var'), vnnlibInfo = []; end
assert(isempty(vnnlibInfo), 'backward-compat fallback to [] failed');

% ----- vnnlibInfo survives a save/load round-trip ------------------------
[~, ~, info] = vnnlib2cora(fullfile(v2Dir, 'multidim_basic.vnnlib'));
save(matFile, 'info');
clear info; S = load(matFile, 'info'); info = S.info;
assert(strcmp(info.version, '2.0'), 'save/load: version lost');
assert(strcmp(info.networks.inputs(1).name, 'X'), 'save/load: input name lost');
% and the reloaded struct still drives the writer correctly
content = aux_write(testDir, (1:6)'/10, [1;2], info);
assert(contains(content, 'X real [2,3]'), 'save/load: reloaded info writes wrong header');

% ----- numeric precision of the %g writer --------------------------------
info = struct('version', '2.0', 'networks', struct( ...
    'inputs', struct('name', 'X', 'dtype', 'real', 'shape', 5), ...
    'output', struct('name', 'Y', 'dtype', 'real', 'shape', 1)));
x_ = [1e-6; 1e6; -3.14159; 0; 1.23456789]; y_ = -1e-10;
content = aux_write(testDir, x_, y_, info);
[xp, ~] = aux_parseV2(content, {'X'}, {'Y'});
assert(numel(xp) == 5, 'precision: not all values written');
assert(abs(xp(1) - 1e-6) <= 1e-12, 'precision: tiny value lost');
assert(abs(xp(2) - 1e6)  <= 1,     'precision: large value lost');
assert(xp(4) == 0, 'precision: zero not preserved');

fprintf('testnn_neuralNetwork_verify_vnncomp_regression: all checks passed.\n');

end


% Auxiliary functions -----------------------------------------------------

function content = aux_write(testDir, x_, y_, vnnlibInfo)
% write a counterexample to a temp file and return its content
tf = fullfile(testDir, 'ce.txt');
fid = fopen(tf, 'w');
writeCounterexample(fid, x_, y_, vnnlibInfo);
fclose(fid);
content = fileread(tf);
end

function [x_, y_] = aux_parseV1(content)
% parse legacy ((X_i v) ... (Y_j v) ...) output
xt = regexp(content, '\(X_\d+\s+([-\d.eE+]+)\)', 'tokens');
yt = regexp(content, '\(Y_\d+\s+([-\d.eE+]+)\)', 'tokens');
x_ = cellfun(@(t) str2double(t{1}), xt)';
y_ = cellfun(@(t) str2double(t{1}), yt)';
end

function [x_, y_] = aux_parseV2(content, inNames, outNames)
% parse v2 output; bucket values by whether the header name is an input or
% output, concatenating in the order they appear (= joint-space order)
lines = strsplit(strtrim(content), newline);
x_ = []; y_ = []; isInput = false; active = false;
for i = 2:numel(lines)             % skip the leading 'sat'
    line = strtrim(lines{i});
    if isempty(line), continue; end
    hdr = regexp(line, '^(\S+)\s+(\w+)\s+\[[\d,]*\]$', 'tokens', 'once');
    if ~isempty(hdr)
        name = hdr{1};
        if any(strcmp(name, inNames))
            isInput = true; active = true;
        elseif any(strcmp(name, outNames))
            isInput = false; active = true;
        else
            active = false;        % unknown variable -> ignore its values
        end
    elseif active
        v = str2double(line);
        if ~isnan(v)
            if isInput, x_ = [x_; v]; else, y_ = [y_; v]; end %#ok<AGROW>
        end
    end
end
end

function names = aux_headerOrder(content)
% return the variable names of header lines in the order they appear
lines = strsplit(strtrim(content), newline);
names = {};
for i = 1:numel(lines)
    hdr = regexp(strtrim(lines{i}), '^(\S+)\s+\w+\s+\[[\d,]*\]$', ...
        'tokens', 'once');
    if ~isempty(hdr), names{end+1} = hdr{1}; end %#ok<AGROW>
end
end

% ------------------------------ END OF CODE ------------------------------
