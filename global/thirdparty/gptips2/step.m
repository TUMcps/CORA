function y = step(x)
%STEP Node. Threshold function that returns 1 if the argument is >= 0 and 0 otherwise.
%
%   Y = STEP(X)
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
%   See also THRESH, IFLTE, MAXX, MINX, NEG, LTH, GTH, GPAND, GPNOT, GPOR

y = ( (x>=0) .* 1) + ( (x<0) .* 0);

