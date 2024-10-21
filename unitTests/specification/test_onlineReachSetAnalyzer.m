function res = test_onlineReachSetAnalyzer
% test_onlineReachSetAnalyzer - unit test function of onlineReachSetAnalyzer
%
% Syntax:
%    res = test_onlineReachSetAnalyzer
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
% See also: -

% Authors:       Florian Lercher
% Written:       21-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% atomic predicates
x = stl('x',2);
ap1 = x(1) < 6;
ap2 = x(2) > 10;

% sets
g = 0.25 * eye(2);
TT = zonotope([5;11],g);
TF = zonotope([5;9],g);
FT = zonotope([7;11],g);
FF = zonotope([7;9],g);
UT = zonotope([6;11],g);
UF = zonotope([6;9],g);
UU = zonotope([6;10],g);
TU = zonotope([5;10],g);
FU = zonotope([7;10],g);

% test propagation frequency
analyzer = onlineReachSetAnalyzer(finally(ap1,stlInterval(0,1)),3);
analyzer.observeSet(TT,stlInterval(0,1,true,false),0);
assert(analyzer.getVerdict() == fourValued.Inconclusive);% observation has not been propagated yet
analyzer.observeSet(FF,stlInterval(1,2,true,false),0);
assert(analyzer.getVerdict() == fourValued.Inconclusive);% observation has not been propagated yet
analyzer.observeSet(UU,stlInterval(2,3,true,false),0);
assert(analyzer.getVerdict() == fourValued.True);% after two more observations the signals are propagated

% test force propagation
analyzer = onlineReachSetAnalyzer(finally(ap1,stlInterval(0,1)),3);
analyzer.observeSet(TT,stlInterval(0,1,true,false),0);
assert(analyzer.getVerdict() == fourValued.Inconclusive);% observation has not been propagated yet
analyzer.forcePropagation();
assert(analyzer.getVerdict() == fourValued.True);% after force propagation the signals are propagated

% test propagation
analyzer = onlineReachSetAnalyzer(until(ap1,globally(ap2,stlInterval(0,2,false)),stlInterval(0,inf)),1);
analyzer.observeSet(TT,stlInterval(0,1,true,false),0);
assert(analyzer.getVerdict() == fourValued.Inconclusive);
analyzer.observeSet(TU,stlInterval(1,2,true,false),0);
assert(analyzer.getVerdict() == fourValued.Inconclusive);
analyzer.observeSet(FT,stlInterval(2,3,true,false),0);
assert(analyzer.getVerdict() == fourValued.Inconclusive);
analyzer.observeSet(UT,stlInterval(3,4,true,false),0);
assert(analyzer.getVerdict() == fourValued.True);

% test observe set
analyzer = onlineReachSetAnalyzer(ap1,1);
analyzer.observeSet(TF,stlInterval(0,1),0);
assert(analyzer.getVerdict() == fourValued.True);
analyzer = onlineReachSetAnalyzer(ap1,1);
analyzer.observeSet(UF,stlInterval(0,1),0);
assert(analyzer.getVerdict() == fourValued.Unknown);
analyzer = onlineReachSetAnalyzer(ap1,1);
analyzer.observeSet(FU,stlInterval(0,1),0);
assert(analyzer.getVerdict() == fourValued.False);

% test plot
analyzer = onlineReachSetAnalyzer(until(ap1,globally(ap2 & ~finally(ap1,stlInterval(3,4,false)),stlInterval(0,inf)),stlInterval(0,1)),3);
figure;

analyzer.plot();
han = analyzer.plot(); %#ok<NASGU>

close;

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
