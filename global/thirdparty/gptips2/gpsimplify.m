function symOut = gpsimplify(symIn,numSteps,verOld,fastMode,seconds)
%GPSIMPLIFY Simplify SYM expressions in a less glitchy way than SIMPLIFY or SIMPLE.
%
%   SYMOUT = GPSIMPLIFY(SYMIN, NUMSTEPS, VEROLD, FAST) simplifies SYMIN
%   using a maximum of NUMSTEPS using the MuPAD engine directly if VEROLD
%   is FALSE and falling back to the standard MATLAB SIMPLIFY method is
%   VEROLD is TRUE.
%
%   The slower (but 'better') MuPAD Engine symbolic:Simplify method is used
%   if FASTMODE is FALSE and the possibly faster symbolic:simplify method
%   is used if FASTMODE = TRUE. The default is FASTMODE = FALSE.
%
%   According to Mathworks, the MuPAD symbolic:simplify method is "faster,
%   but less reliable and controllable" than the algorithm used by MuPAD
%   symbolic:Simplify.
%
%   It is recommended that you use the functions GPPRETTY and GPMODEL2SYM
%   to simplify GPTIPS expressions, rather than calling this function
%   directly.
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also GPPRETTY, GPMODEL2SYM, GPMODEL2STRUCT, SYM/VPA

if nargin < 4 || isempty(fastMode)
    fastMode = false;
end

if nargin < 5 || isempty(seconds)
    seconds = '0.1';
end

try
    if verOld %for 2007 MATLAB, default to standard simplify method
        symOut = simplify(symIn);
    else
        %otherwise call Mupad symengine directly
        
        if fastMode
            %Note that according to:
            %http://www.mathworks.co.uk/help/symbolic/mupad_ug/use-general-simplification-functions.html
            %The MuPAD symbolic:simplify method "searches for a simpler
            %form by rewriting the terms of an expression." and "generally,
            %this method is faster, but less reliable and controllable than
            %the algorithm used by symbolic:Simplify."
            symOut = evalin(symengine,['symbolic:simplify(' char(symIn) ')']);
        else %and similarly, symbolic:Simplify "performs more extensive search
            %for a simple form of an expression. For some expressions, this
            %function can be slower, but more flexible and more powerful
            %than simplify."
            symOut = evalin(symengine,['symbolic:Simplify(' char(symIn) ',Seconds = ' num2str(seconds) ', Steps = ' num2str(numSteps) ')']);
        end
    end
catch
    symOut = symIn;
end