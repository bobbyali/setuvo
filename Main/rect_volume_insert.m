% =========================================================================
% rect_volume_insert.m
% Rehan Ali, 9th December 2011
%
% Insert subvolume into bigger volume using rect coords from previous 
% user selection.
%
% Input:
%   I    -    input subvolume
%   rect -    coordinates selected by user
%   z    -    optional extra z frames (defined as vector [a:b])
%   
% Output:
%   J    -    final volume
%
% =========================================================================  

function [J] = rect_volume_insert(J,I,rect,z_start,z_end)

    rect = floor(rect);
    if nargin < 4
        J(rect(2):rect(2)+rect(4), rect(1):rect(1)+rect(3),:) = I;
    else
        J(rect(2):rect(2)+rect(4), rect(1):rect(1)+rect(3),z_start:z_end) = I;
    end
    