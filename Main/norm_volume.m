% =========================================================================
% norm_volume.m
% Rehan Ali, 8th December 2011
%
% Normalises volume intensity values to 0-1.
% =========================================================================

function W = norm_volume(V)
    
    W = V - min(min(min(V)));
    W = W / max(max(max(W)));