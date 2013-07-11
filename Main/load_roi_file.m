% =========================================================================
% load_roi_file.m
% Rehan Ali, 23rd May 2012
%
% Load an RT_image ROI file and crop it down to a particular size to
% compare with tumor segmentation matrices.
% =========================================================================

function roi_vol_tumor = load_roi_file(filePath,rect_tumor,startSlice,endSlice,dx,dy,dz)
  
    if filePath(end-2:end) == 'mat'
        % mat files from Amira
        load(filePath);
        slashes  = strfind(filePath,'/');
        var_name = filePath(slashes(end)+1:end-11);
        var_name = strrep(var_name,'.','_');
        filename = ['Amira_' var_name '_labels_mat']; 
        eval(['roi_vol = ' filename ';']);
        eval(['clear ' filename]);
        roi_vol = permute(roi_vol,[3 1 2]);
        roi_vol = flipdim(roi_vol,1);
        roi_vol = flipdim(roi_vol,3);
    else
        % bin/raw files from RT_image
        fh   = fopen(filePath,'r','l');
        bin  = fread(fh,'uchar');
        fclose(fh);
        bin(1:8) = [];
        roi_vol  = reshape(bin,[dy dx dz]);
        clear bin
        for i = 1 : dz
            roi_vol1(:,:,i) = rotateCT(roi_vol(:,:,i));
        end
        clear roi_vol;
        roi_vol = roi_vol1;
        clear roi_vol1;
        roi_vol = uint8(roi_vol);
    end
    roi_vol_tumor = rect_volume_select(roi_vol,rect_tumor,startSlice,endSlice);
    roi_vol_tumor(roi_vol_tumor > 0) = 1;   
    clear roi_vol
    
    