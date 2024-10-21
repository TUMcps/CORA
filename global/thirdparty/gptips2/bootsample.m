function xind = bootsample(x,sampleSize)
%BOOTSAMPLE Get an index vector to sample with replacement a data matrix X.
%
%   XIND = BOOTSAMPLE(X,SAMPLESIZE) samples the data matrix X to return the
%   row indices XIND of X.
% 
%   Copyright (c) 2009-2015 Dominic Searson 
%
%   GPTIPS 2

xSize = size(x,1);

%generate the required sample indices and sample
xind = ceil(xSize*rand(sampleSize,1));


