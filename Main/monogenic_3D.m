% EXTENSION TO 3D

% =========================================================================
% monogenic.m
% Dr Matthew Mellor and Weiwei Zhang, University of Oxford
% Modified by Rehan Ali, 19th October 2007
% Modified by Tuende Szilagyi, July, 2008 (FT sampling corrected)
% Modified by Tuende Szilagyi, January, 2009 (Speed up- allocate
% memory,Instantaneous frequency added)
% Modified by Rehan Ali, 23 May 2012, fixed bug that caused code to crash with
% non-square input images (see lines appended with %! symbols)
%
% Generates a monogenic signal representation of a 2D image, which is used
% to compute local energy, local phase, and local orientation. It uses a
% specific scale-invariant filter (which can be defined as bandpass for
% feature detection, or lowpass with zero DC for approximating an inverse
% laplacian). Local amplitude can be seen as a descriptor of the strength
% of a local feature, whilst local phase describes the type of feature
% present. Local orientation is a measure which points in the direction of
% maximal signal variance. 
%
%
% The monogenic signal framework and theory was initially developed by
% Felsberg and Sommer:
%
% "The monogenic signal", M Felsberg and G Sommer, IEEE
% Transactions on Signal Processing, Dec 2001, Volume 49, Issue 12,
% 3136-3144
%
%
% For details of the Mellor-Brady filter and its application, refer to the
% following references:
%
% "Phase mutual information as a similarity measure for registration",
% Matthew Mellor and Michael Brady, Medical Image Analysis, 
% Volume 9, Issue 4, August 2005, Pages 330-343 
%
% "Phase Methods for Non-Rigid Medical Image Registration", Matthew Mellor,
% PhD Thesis, University of Oxford, 2004
%
%
% For details on how to use this filter to perform an inverse laplacian,
% see the reference:
%
% "Advanced Phase-Based Segmenation of Multiple Cells with Brightfield
% Microscopy",   R. Ali, M. Gooding, M. Christlieb, J.M. Brady, in Proc.
% International Symposium   on Biomedical Imaging, May 2008
%
%
% Input:
%   image    -    input image
%   alpha    -    alpha parameter for Mellor-Brady filter (optional)
%   beta     -    beta parameter for Mellor-Brady filter  (optional)
%
% Output:
%   f        -    bandpassed image
%   le       -    local energy image
%   lo       -    local orientation image
%   lp       -    local phase image
%
%
% =========================================================================
function  [monogenicStack, le, lp] = monogenic_3D(image,NN, alpha,beta)

    warning off

    [dimx, dimy, dimz] = size(image);
  
    
    if nargin == 2

      % set up standard Mellor-Brady filter params for feature detection if
      % none are specifiied
      alpha = 2.25;
      beta  = 0.25; 
 
    end;
    
    % set up zero-centered grid to compute values on
    maxu = max(dimx, dimy); maxu = max(maxu, dimz); 
    % where dimx and dimy are the input image coordinates
    %! [u1, u2, u3] = meshgrid(linspace(-maxu, maxu, dimx*2), linspace(-maxu, maxu, dimy*2), linspace(-maxu, maxu, dimz*2));
    [u1, u2, u3] = meshgrid(linspace(-maxu, maxu, dimy*2), linspace(-maxu, maxu, dimx*2), linspace(-maxu, maxu, dimz*2));

    hyp = sqrt(u1.^2 + u2.^2 + u3.^2);

    
    for N = 1: NN
        
    % set up Riesz Transform filters in 3D
    %! v1 = ([0:dimx-1]-round(dimx/2)); v2 = ([0:dimy-1]-round(dimy/2)); v3 = ([0:dimz-1]-round(dimz/2));
    v1 = ([0:dimy-1]-round(dimy/2)); v2 = ([0:dimx-1]-round(dimx/2)); v3 = ([0:dimz-1]-round(dimz/2));
    [v1 v2 v3] = meshgrid(v1, v2, v3);
    
    % set up the Riesz filters
    hyp2 = sqrt(v1.^2 + v2.^2 + v3.^2);
    V1 = i*v1./(hyp2 + eps); V2 = i*v2./(hyp2 + eps); V3 = i*v3./(hyp2 + eps); 
       
    % upsample image in 3D no built-in Matlab function for this
    d = size(image);
    %! [xi, yi, zi] = meshgrid(linspace(1, d(1), 2*d(1)), linspace(1, d(2), 2*d(2)), linspace(1, d(3), 2*d(3)));
    [xi, yi, zi] = meshgrid(linspace(1, d(2), 2*d(2)), linspace(1, d(1), 2*d(1)), linspace(1, d(3), 2*d(3)));
    img_resized = interp3(image, xi, yi, zi, 'linear'); % you could change the interpolation method here
        
    fftimg = fftn(img_resized);

    % generate Mellor-Brady filter 3D (kernels in domain)
    kern_p = hyp.^-(alpha + beta*(N+1)); 
    kern_p = kern_p/sum(sum(sum(kern_p)));
    kern_n = hyp.^-(alpha + beta*N);
    kern_n = kern_n/sum(sum(sum(kern_n)));
    kern = kern_p - kern_n;
    

    % convolve image with filter
    im_r = real(ifftn(fftimg.*fftn(ifftshift(kern))));

    % subsample image back to original size
    for x = 1:dimx
        for y = 1:dimy
            for z = 1:dimz
                %imf(x, y, z) = im_r(x*2-1, y*2-1, z*2-1);
                imf(x, y, z) = im_r(x*2-1, y*2-1, z*2-1);
            end    
        end
    end
  
    % generate spherical quadrature filter triplet
    f(:, :, :, N) = imf;
    monogenicStack(:, :, :, 1, N) = imf; 
      
    fftimg3 = fftshift(fftn(imf)); 
    monogenicStack(:, :, :, 2, N) = real(ifftn(ifftshift(V1.*fftimg3)));
    monogenicStack(:, :, :, 3, N) = real(ifftn(ifftshift(V2.*fftimg3)));
    monogenicStack(:, :, :, 4, N) = real(ifftn(ifftshift(V3.*fftimg3)));

    % generate band-passed image
    f(:, :, :, N) = imf;                                                              
                         
    % generate local energy image
    le(:, :, :, N)= sqrt(monogenicStack(:, :, :, 1,N).^2 + monogenicStack(:, :, :, 2, N).^2 +...
               monogenicStack(:, :, :, 3, N).^2 + monogenicStack(:, :, :, 4, N).^2);

    % generate local phase image 3D
    lp(:, :, :, N)= atan2(monogenicStack(:, :, :, 1, N), sqrt(monogenicStack(:, :, :, 2, N).^2 + monogenicStack(:, :, :,3, N).^2 + monogenicStack(:, :, :, 4, N).^2)); 

end