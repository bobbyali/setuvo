% =========================================================================
% downsize3D.m
% Rehan Ali, 10th March 2012
%
% Simple downsampling algorithm that halves the number of row, columns and
% planes in a 3D volume.
% =========================================================================
function W = downsize3D(V)

    [dimX dimY dimZ] = size(V);
    
    newX = floor(dimX/2);
    newY = floor(dimY/2);
    newZ = floor(dimZ/2);
    
    W = zeros(newX,newY,newZ);
    
    for i = 1:newX
        for j = 1:newY
            for k = 1:newZ
                W(i,j,k) = V(2*i,2*j,2*k);                               
            end
        end
    end