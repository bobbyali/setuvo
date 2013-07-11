% =========================================================================
% upsize3D.m
% Rehan Ali, 10th March 2012
%
% Simple upsampling algorithm that repeats the rows, columns and planes in
% a 3D volume. (Nearest neighbour interpolation) Only to be used for
% integer maps, as it'll pixelise any image data.
% =========================================================================
function C = upsize3D(V,newX,newY,newZ)

    [dimX dimY dimZ] = size(V);
       
    A = zeros(dimX,newY,dimZ);
    B = zeros(newX,newY,dimZ);
    C = zeros(newX,newY,newZ);
    
    j = 1;
    for i = 1:dimY
        A(1:dimX,j,1:dimZ) = V(:,i,:);
        j = j + 1;
        A(1:dimX,j,1:dimZ) = V(:,i,:);
        j = j + 1;
    end
    
    j = 1;
    for i = 1:dimX
        B(j,1:newY,1:dimZ) = A(i,:,:);
        j = j + 1;
        B(j,1:newY,1:dimZ) = A(i,:,:);
        j = j + 1;
    end
    
    j = 1;
    for i = 1:dimZ
        C(:,:,j) = B(:,:,i);
        j = j + 1;
        C(:,:,j) = B(:,:,i);
        j = j + 1;
    end
    