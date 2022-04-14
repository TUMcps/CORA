function [ corapath ] = coraroot()
%CORAROOT
%   Returns the CORA root path

s = which('coraroot');
corapath = fileparts(fileparts(fileparts(s)));