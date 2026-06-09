function res = test_nn_converter_vnnlib2cora_equivalence()
% test_nn_converter_vnnlib2cora_equivalence - parser-version equivalence:
%    for each canonical 1.x fixture, auto-convert to 2.0 in memory, parse
%    both versions through the dispatcher, and assert that X0 and spec are
%    identical. The 2.0 parser is correct for legacy specs iff this passes.
%
% Syntax:
%    res = test_nn_converter_vnnlib2cora_equivalence()
%
% Inputs:
%    -
%
% Outputs:
%    res - true on success
%
% See also: vnnlib2cora

% Authors:       Benedikt Kellner
% Written:       19-April-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

baseDir = [CORAROOT '/models/Cora/nn/unitTests/vnnlib'];
files = { ...
    'axas_xu_prop_3.vnnlib', ...
    'mnistfc_prop_11_0.05.vnnlib', ...
    'rl_benchmark_dubinsrejoin_case_unsafe_11.vnnlib'};

for i = 1:numel(files)
    fname = fullfile(baseDir, files{i});
    aux_assertEquivalent(fname);
end

% Also exercise a few 1.x ACAS-Xu props from models/Cora/nn/
extraDir = [CORAROOT '/models/Cora/nn'];
extraFiles = {'prop_1.vnnlib','prop_2.vnnlib','prop_5.vnnlib', ...
              'prop_7.vnnlib','prop_8.vnnlib','prop_9.vnnlib'};
for i = 1:numel(extraFiles)
    fname = fullfile(extraDir, extraFiles{i});
    if isfile(fname)
        aux_assertEquivalent(fname);
    end
end

res = true;

end


% Auxiliary functions -----------------------------------------------------

function aux_assertEquivalent(fname)
[X0_v1, spec_v1] = vnnlib2cora(fname);

tmpFile = [tempname, '.vnnlib'];
cleaner = onCleanup(@() aux_safeDelete(tmpFile));
fid = fopen(tmpFile, 'w');
fwrite(fid, aux_convertVnnlib1To2(fileread(fname)));
fclose(fid);

[X0_v2, spec_v2] = vnnlib2cora(tmpFile);
clear cleaner; %#ok<CLEAR>

[eq, reason] = aux_vnnlibResultsEqual(X0_v1, spec_v1, X0_v2, spec_v2);
assert(eq, '%s: %s', fname, reason);

end

function aux_safeDelete(tmpFile)
if isfile(tmpFile)
    delete(tmpFile);
end
end

function text2 = aux_convertVnnlib1To2(text1)
% In-memory converter from VNN-LIB 1.x to 2.0 spec text.
% Mirrors the official to_vnnlib2.py logic for single-network specs:
%   count max X_<n> and Y_<n> declare-const indices, strip declare-const
%   lines, prepend (vnnlib-version <2.0>) header with a synthetic
%   declare-network block, rewrite X_<n> -> X[<n>] and Y_<n> -> Y[<n>].
xToks = regexp(text1, '\(\s*declare-const\s+X_(\d+)\s+Real\s*\)', 'tokens');
yToks = regexp(text1, '\(\s*declare-const\s+Y_(\d+)\s+Real\s*\)', 'tokens');

nX = 0;
for i = 1:numel(xToks), nX = max(nX, str2double(xToks{i}{1}) + 1); end
nY = 0;
for i = 1:numel(yToks), nY = max(nY, str2double(yToks{i}{1}) + 1); end

if nX == 0 || nY == 0
    throw(CORAerror('CORA:converterIssue', ...
        'convertVnnlib1To2: no declare-const X or Y lines found'));
end

body = regexprep(text1, '\(\s*declare-const\s+\w+\s+Real\s*\)\s*\r?\n?', '');
body = regexprep(body, '\<X_(\d+)\>', 'X[$1]');
body = regexprep(body, '\<Y_(\d+)\>', 'Y[$1]');

header = sprintf(['(vnnlib-version <2.0>)\n\n' ...
                  '(declare-network N\n' ...
                  '    (declare-input  X real [%d])\n' ...
                  '    (declare-output Y real [%d])\n' ...
                  ')\n\n'], nX, nY);
text2 = [header, body];
end

function [eq, reason] = aux_vnnlibResultsEqual(X0a, specA, X0b, specB)
% Structural equality check on a (X0, spec) pair returned by vnnlib2cora.
eq = false;

if numel(X0a) ~= numel(X0b)
    reason = sprintf('X0 cell size mismatch (%d vs %d)', numel(X0a), numel(X0b));
    return;
end
for i = 1:numel(X0a)
    if ~strcmp(class(X0a{i}), class(X0b{i}))
        reason = sprintf('X0{%d} class mismatch (%s vs %s)', i, ...
            class(X0a{i}), class(X0b{i}));
        return;
    end
    if isa(X0a{i}, 'interval')
        if ~isequal(infimum(X0a{i}), infimum(X0b{i})) ...
                || ~isequal(supremum(X0a{i}), supremum(X0b{i}))
            reason = sprintf('X0{%d} bounds mismatch', i);
            return;
        end
    end
end

if length(specA) ~= length(specB)
    reason = sprintf('spec length mismatch (%d vs %d)', length(specA), length(specB));
    return;
end
for i = 1:length(specA)
    if ~strcmp(specA(i).type, specB(i).type)
        reason = sprintf('spec(%d).type mismatch (%s vs %s)', i, ...
            specA(i).type, specB(i).type);
        return;
    end
    if ~isequal(specA(i).set.A, specB(i).set.A)
        reason = sprintf('spec(%d).set.A mismatch', i);
        return;
    end
    if ~isequal(specA(i).set.b, specB(i).set.b)
        reason = sprintf('spec(%d).set.b mismatch', i);
        return;
    end
end

eq = true;
reason = '';
end

% ------------------------------ END OF CODE ------------------------------
