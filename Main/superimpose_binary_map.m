% =========================================================================
% superimpose_binary_map.m
% Rehan Ali, 19th October 2009
%
% Superimposes 1-3 binary maps over grayscale image.
%
% Inputs:
%       I           Input image, m x n, scaled to 0-1
%       B1          Binary map, m x n
%       B2          Binary map, m x n
%       B3          Binary map, m x n
%
% Outputs:
%       J           Final rgb image
%
% v1.1  Modified from original which required user to generate m x n x 3
% matrix of binaries to superimpose
% =========================================================================

function J = superimpose_binary_map(I,B1,B2,B3)

    if nargin == 2
        B2 = [];
        B3 = [];
    elseif nargin == 3
        B3 = [];
    end

    [dimX dimY] = size(I);
    J = zeros(dimX,dimY,3);
    
    % normalise input image
    if max(max(I)) > 1
        I = I / max(max(I));
    end
    
    % binarise image maps if not already binary
    if ~isempty(B1) && max(max(B1)) > 1
        B1(B1 > 0) = 1;
    end
    
    if ~isempty(B2) && max(max(B2)) > 1
        B2(B2 > 0) = 1;
    end
    
    if ~isempty(B3) && max(max(B3)) > 1
        B3(B3 > 0) = 1;
    end
    
    J(:,:,1) = I;
    J(:,:,2) = I;
    J(:,:,3) = I;
       
    % apply binaries to image for red, green and blue channels
    if ~isempty(B1)
        temp = I;
        temp(B1 == 1) = 1;
        J(:,:,1) = temp;
    end
    
    if ~isempty(B2)
        temp = I;
        temp(B2 == 1) = 1;
        J(:,:,2) = temp;
    end
    
    if ~isempty(B3)
        temp = I;
        temp(B3 == 1) = 1;
        J(:,:,3) = temp;
    end
    
    clear I B1 B2 B3 temp
    