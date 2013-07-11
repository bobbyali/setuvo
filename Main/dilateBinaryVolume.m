% =========================================================================
% dilateBinaryVolume.m
% Rehan Ali, 9th December 2011
%
% Takes a binary volume mask, where bg = 0 and obj = 1, and expands each 
% object region by one voxel in each direction.
%
% Currently doesn't treat object voxels on boundary.
% =========================================================================

function J = dilateBinaryVolume(I)

    [dx dy dz] = size(I);
    J = I;
    for i = 2 : dx-1
        for j = 2 : dy-1
            for k = 2 : dz-1
                if I(i,j,k) == 1
                    J(i-1,j,k) = 1;
                    J(i+1,j,k) = 1;
                    J(i,j-1,k) = 1;
                    J(i,j+1,k) = 1;
                    J(i,j,k-1) = 1;
                    J(i,j,k+1) = 1;
                    
                    J(i-1,j-1,k-1) = 1;
                    J(i-1,j  ,k-1) = 1;
                    J(i-1,j+1,k-1) = 1;
                    
                    J(i+1,j-1,k-1) = 1;
                    J(i+1,j  ,k-1) = 1;
                    J(i+1,j+1,k-1) = 1;
                    
                    J(i-1,j-1,k+1) = 1;
                    J(i-1,j  ,k+1) = 1;
                    J(i-1,j+1,k+1) = 1;
                    
                    J(i+1,j-1,k+1) = 1;
                    J(i+1,j  ,k+1) = 1;
                    J(i+1,j+1,k+1) = 1;
                end
            end
        end
    end
                    