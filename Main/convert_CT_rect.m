% =========================================================================
% convert_CT_rect.m
% Rehan Ali, 22nd May 2012
%
% Takes a rectangle ROI drawn using getrect on a 90 degree rotated image
% (produced by the rotateCT.m function), and rotates it back to map onto
% the original CT image. Need the dx dimension (size(X,1)) for the
% original unrotated image.
% =========================================================================

function newrect = convert_CT_rect(rect,dx)

    newrect = [rect(2) dx-rect(1)-rect(3) rect(4) rect(3)];