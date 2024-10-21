function y = neg(x)
%NEG Node. Returns -1 times the argument.
%
%   Y = NEG(X)
%
%   Remarks:
%
%   Computed on an element by element basis, and if X is scalar then Y will
%   be scalar.
%
%   Copyright (c) 2009-2015 Dominic Searson 
%   
%   GPTIPS 2
%
%   See also GPNOT, GPAND, GPOR, MAXX, MINX, STEP, THRESH, IFLTE

y = x .* -1;

