% =========================================================================
% dice.m
% Rehan Ali, 22nd May 2012
%
% Computes the DICE coefficient between two binary masks (2D or 3D)
% representing two segmentation results. Typically used to compare an
% automatic segmentation vs a manual segmentation ground truth.
% 
% http://sve.loni.ucla.edu/instructions/metrics/dice/?wscr=1680x1050
% =========================================================================

function D = dice(seg1,seg2)

    overlap = seg1 + seg2;
    AuB = sum(sum(sum(overlap == 2)));
    A   = sum(sum(sum(seg1)));
    B   = sum(sum(sum(seg2)));
    %disp(['Vol A = ' num2str(A) ', B = ' num2str(B)]);
    D = 2*AuB / (A + B);