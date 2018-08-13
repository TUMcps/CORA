function [ corapath ] = coraroot()
%CORAROOT
%   Returns the CORA root path

s = which('coraroot');
% s is something, such as C:\Repositories\Matlab\CORA\global functions\coraroot.m
% remove the last two parts
c = regexp(s, '[\\/]');
if (size(c) < 2)
    error('path not found');
end

corapath = s(1:c(end-1)-1);


% % Guesses the cora root path by the include path
% corapath = '';
% for i=length(p):-1:1
%     if length(strfind(p{i}, 'CORA')) > 0 && (length(corapath) == 0 || length(corapath) > length(p{i}))
%         corapath = p{i};
%     end
%     
% end
% 
% if length(corapath) == 0
%     error('CORA path not found!');
% end
