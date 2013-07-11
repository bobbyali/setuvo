% =========================================================================
% rect_volume_select.m
% Rehan Ali, 9th December 2011
%
% Select part of a volume based on a previous user manual selection, using
% the coordinates saved from the user operation
%
% Input:
%   I    -    input image volume
%   rect -    coordinates selected by user
%   z    -    optional extra start-end z frames (as vector [a:b])
%   
% Output:
%   J    -    final subvolume
%
% =========================================================================  

function [J] = rect_volume_select(I,rect,z_start,z_end)

    rect = floor(rect);

    if nargin < 3
        J = I(rect(2):rect(2)+rect(4), rect(1):rect(1)+rect(3),:);     
    else
        J = I(rect(2):rect(2)+rect(4), rect(1):rect(1)+rect(3),z_start:z_end);
    end
    