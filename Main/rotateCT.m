% =========================================================================
% rotateCT.m
% Rehan Ali, 9th December 2011
%
% Rotates a slice from a CT image so that it appears the correct way up
% when opened in a matlab figure.
%
% =========================================================================
function J = rotateCT(I)

    J = fliplr(I');