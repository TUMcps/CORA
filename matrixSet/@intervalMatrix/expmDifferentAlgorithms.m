function result = expmDifferentAlgorithms(intMat,maxOrder,styleOfCalculation,varargin)
% returns the approximation of e^intMat using different algorithms 
%    with maxOrder iterations. It is used as a fassade to access the
%    algorithms in the private directory
%
% Syntax:  
%    result = expmDifferentAlgorithms(intMat,maxOrder,styleOfCalculation,varargin);
%
% Inputs:
%    intMat - intervalMatrix (nxn)
%    maxOrder - maximum order of the TaylorSeries, has to be > abs(intMat) +2
%    styleOfCalculation - a number used to indicate the algorithm which
%               should be used for the caluclation of the approximation 
%    varargin - is a list which contains optional parameter for other
%               calculations with extra knowledge needed.
%               e.g. the potential for the scaling and squaring process
% 
% Outputs:
%    result - the exponentiation with the chosen algorithm
%
% Example: 
%
% Other m-files required: hornerTaylorSeries.m, taylorSeries.m,
%           intervalMatrixRemainder.m, scalingSquaringHornerTaylorSeries
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Ivan Brkan
% Written:      23-April-2019
% Last update:  29-April-2019
% Last revision:---

%------------- BEGIN CODE --------------

try

    switch styleOfCalculation
   
        case 0 
            % case of truncated Taylor series
            result = taylorSeries(intMat,maxOrder);
        case 1 
            % case of truncated Taylor series using horner scheme
            result = hornerTaylorSeries(intMat,maxOrder);
        case 2 
            % case of truncated Taylor series using horner scheme and the
            switch length(varargin)
                case 0       
                    % no potential was given as attribute,
                    % so 1 is used as the default value
                    result = scalingSquaringHornerTaylorSeries(intMat,maxOrder,1);
                case 1
                    tmp = cell2mat(varargin);
                    result = scalingSquaringHornerTaylorSeries(intMat,maxOrder,tmp(1));
                otherwise 
                    result = [];
            end

        case 4 
            r = 1;
            result = expm(intMat,r,maxOrder);

        otherwise
            result = [];        
    end
    
catch
    result = [];
end

end

