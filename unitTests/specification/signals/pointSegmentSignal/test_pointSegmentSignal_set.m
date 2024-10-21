function res = test_pointSegmentSignal_set
% test_pointSegmentSignal_set - unit test function of set
%
% Syntax:
%    res = test_pointSegmentSignal_set
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Florian Lercher
% Written:       16-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% signals
tt = true;
ff = false;
zeroT = pointSegmentSignal(0,[tt,tt]);
zeroF = ~zeroT;
lb = 10;
ub = 15;
ccT = pointSegmentSignal([0,lb,ub],[ff,ff,tt,tt,tt,ff]);
ccF = ~ccT;
coT = pointSegmentSignal([0,lb,ub],[ff,ff,tt,tt,ff,ff]);
coF = ~coT;
ocT = pointSegmentSignal([0,lb,ub],[ff,ff,ff,tt,tt,ff]);
ocF = ~ocT;
ooT = pointSegmentSignal([0,lb,ub],[ff,ff,ff,tt,ff,ff]);
ooF = ~ooT;
sigs = {ccT,coT,ocT,ooT};

% test case definition
test_cases = {
    % {sig, int, val, expected}
    % setting blank true signal
    {zeroT, stlInterval(lb,ub,true,true), ff, ccF};
    {zeroT, stlInterval(lb,ub,true,false), ff, coF};
    {zeroT, stlInterval(lb,ub,false,true), ff, ocF};
    {zeroT, stlInterval(lb,ub,false,false), ff, ooF};
    {zeroT, stlInterval(lb,ub,true,true), tt, zeroT};

    % setting blank false signal
    {zeroF, stlInterval(lb,ub,true,true), tt, ccT};
    {zeroF, stlInterval(lb,ub,true,false), tt, coT};
    {zeroF, stlInterval(lb,ub,false,true), tt, ocT};
    {zeroF, stlInterval(lb,ub,false,false), tt, ooT};
    {zeroF, stlInterval(lb,ub,true,true), ff, zeroF};

    % setting at zero
    {zeroT, stlInterval(0), ff, pointSegmentSignal(0,[ff,tt])};
    {zeroT, stlInterval(0,1,true,true), ff, pointSegmentSignal([0,1],[ff,ff,ff,tt])};
    {zeroT, stlInterval(0,1,true,false), ff, pointSegmentSignal([0,1],[ff,ff,tt,tt])};
    {zeroT, stlInterval(0,1,false,true), ff, pointSegmentSignal([0,1],[tt,ff,ff,tt])};
    {zeroT, stlInterval(0,1,false,false), ff, pointSegmentSignal([0,1],[tt,ff,tt,tt])};

    % setting exact interval
    {ccT, stlInterval(lb,ub,true,true), tt, ccT};
    {ccT, stlInterval(lb,ub,true,true), ff, zeroF};
    {coT, stlInterval(lb,ub,true,false), tt, coT};
    {coT, stlInterval(lb,ub,true,false), ff, zeroF};
    {ocT, stlInterval(lb,ub,false,true), tt, ocT};
    {ocT, stlInterval(lb,ub,false,true), ff, zeroF};
    {ooT, stlInterval(lb,ub,false,false), tt, ooT};
    {ooT, stlInterval(lb,ub,false,false), ff, zeroF};

    % setting boundary point left
    {ccT, stlInterval(lb), ff, pointSegmentSignal([0,lb,ub],[ff,ff,ff,tt,tt,ff])};
    {ooT, stlInterval(lb), ff, ooT};
    {ccT, stlInterval(lb,lb+1,true,true), ff, aux_makeOverLeft(ccT,true,true)};
    {ooT, stlInterval(lb,lb+1,true,true), ff, aux_makeOverLeft(ooT,true,true)};
    {ccT, stlInterval(lb,lb+1,false,true), ff, pointSegmentSignal([0,lb,lb+1,ub],[ff,ff,tt,ff,ff,tt,tt,ff])};
    {ooT, stlInterval(lb,lb+1,false,true), ff, aux_makeOverLeft(ooT,false,true)};

    % setting boundary point right
    {ccT, stlInterval(ub), ff, pointSegmentSignal([0,lb,ub],[ff,ff,tt,tt,ff,ff])};
    {ooT, stlInterval(ub), ff, ooT};
    {ccT, stlInterval(ub-1,ub,true,true), ff, aux_makeOverRight(ccT,true,true)};
    {ooT, stlInterval(ub-1,ub,true,true), ff, aux_makeOverRight(ooT,true,true)};
    {ccT, stlInterval(ub-1,ub,true,false), ff, pointSegmentSignal([0,lb,ub-1,ub],[ff,ff,tt,tt,ff,ff,tt,ff])};
    {ooT, stlInterval(ub-1,ub,true,false), ff, aux_makeOverRight(ooT,true,false)};

    % empty interval
    {zeroT, stlInterval(), ff, zeroT};
    {zeroT, stlInterval(), tt, zeroT};
    {zeroF, stlInterval(), ff, zeroF};
    {zeroF, stlInterval(), tt, zeroF};

    % fullspace interval
    {zeroT, stlInterval(0,inf), ff, zeroF};
    {zeroT, stlInterval(0,inf), tt, zeroT};
    {zeroF, stlInterval(0,inf), ff, zeroF};
    {zeroF, stlInterval(0,inf), tt, zeroT};
    {ccT, stlInterval(0,inf), ff, zeroF};
    {coT, stlInterval(0,inf), ff, zeroF};
    {ocT, stlInterval(0,inf), ff, zeroF};
    {ooT, stlInterval(0,inf), ff, zeroF};
};

% setting inside interval
test_cases = [test_cases; aux_genTests(lb+1,ub-1,@aux_makeInner)];

% setting inner point
p = (lb+ub)/2;
test_cases = [test_cases; aux_genTests(p,p,@aux_makeInnerPoint)];

% setting over left boundary
test_cases = [test_cases; aux_genTests(lb-1,lb+1,@aux_makeOverLeft)];

% setting over right boundary
test_cases = [test_cases; aux_genTests(ub-1,ub+1,@aux_makeOverRight)];

% run tests
for i = 1:length(test_cases)
    sig = test_cases{i}{1};
    int = test_cases{i}{2};
    val = test_cases{i}{3};
    expected = test_cases{i}{4};
    actual = sig.set(int,val);
    assertLoop(actual == expected,i)
end

res = true;


% Auxiliary functions -----------------------------------------------------

function tests = aux_genTests(lower,upper,oracle)
    tests = cell(length(sigs)*4,1);
    for j = 1:length(sigs)
        s = sigs{j};
        jj = (j-1)*4+1;
        tests{jj} = {s, stlInterval(lower,upper,true,true), ff, oracle(s,true,true)};
        tests{jj+1} = {s, stlInterval(lower,upper,true,false), ff, oracle(s,true,false)};
        tests{jj+2} = {s, stlInterval(lower,upper,false,true), ff, oracle(s,false,true)};
        tests{jj+3} = {s, stlInterval(lower,upper,false,false), ff, oracle(s,false,false)};
    end
end

function inner = aux_makeInner(sig,lc,rc)
    tp = sig.timePoints;
    val = sig.values;
    newTp = [lb+1,ub-1];
    newVal = [~lc,ff,~rc,tt];
    inner = pointSegmentSignal([tp(1:2),newTp,tp(3:end)],[val(1:4),newVal,val(5:end)]);
end

function inner = aux_makeInnerPoint(sig,lc,rc)
    tp = sig.timePoints;
    val = sig.values;
    if lc && rc
        newTp = p;
        newVal = [ff,tt];
        inner = pointSegmentSignal([tp(1:2),newTp,tp(3:end)],[val(1:4),newVal,val(5:end)]);
    else
        inner = sig;
    end
end

function inner = aux_makeOverLeft(sig,~,rc)
    tp = sig.timePoints;
    val = sig.values;
    newTp = lb + 1;
    newVal = [~rc,tt];
    inner = pointSegmentSignal([tp(1:1),newTp,tp(3:end)],[val(1:2),newVal,val(5:end)]);
end

function inner = aux_makeOverRight(sig,lc,~)
    tp = sig.timePoints;
    val = sig.values;
    newTp = ub - 1;
    newVal = [~lc,ff];
    inner = pointSegmentSignal([tp(1:2),newTp,tp(4:end)],[val(1:4),newVal,val(7:end)]);
end

end

% ------------------------------ END OF CODE ------------------------------
