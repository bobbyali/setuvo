% =========================================================================
% segment_batch_item.m
% Rehan Ali, 23rd May 2012
%
% Executes a batched segmentation task for CT tumor segmentation.
% =========================================================================

function segment_batch_item(itemPath)

    % load seg task data
    load(itemPath);
    seg = uint8(zeros(size(app.vol.lp))); % saves seg results in viewable orientation    
    
    % run level set
    tic
    [seg,phi,ls_vols,tmap] = levelset3D(app.vol.lp, app.vol.init, app.maxits, app.eta, app.T, app.alpha);    
    tSeg = toc;
    
    % start post-processing
    seg  = resize(seg,[app.fulldx app.fulldy app.fulldz]);
    phi  = resize(phi,[app.fulldx app.fulldy app.fulldz]);
    tmap = resize(tmap,[app.fulldx app.fulldy app.fulldz]);
    
    app.vol.seg = seg;
    app.vol.phi = phi;
    app.vol.tmap = tmap;
    app.ls_vols = ls_vols;
    app.stage = 2;    
    
    i0 = max(strfind(itemPath,'/'));
    i1 = length(itemPath)-4;
    strResultsLabel = itemPath(i0:i1);
    save(['results/' strResultsLabel '.mat'],'app');
    
    % do plot of seg volume vs iterations
    figure;
    plot(app.ls_vols);
    xlabel('Iterations');ylabel('Volume');title('Size of Zero Level Set over Time');
    saveas(gcf,['results/' strResultsLabel '_convergence_plot.png'],'png');

    % create movie flythrough of 3D segmentation
    figure;
    vidObj = VideoWriter(['results/' strResultsLabel '_3Dflythru.avi']);
    open(vidObj);
    for i = 1 : app.dz

        fullimg = superimpose_binary_map(rotateCT(app.vol.tumor(:,:,i)),rotateCT(app.vol.seg(:,:,i)),rotateCT(app.vol.init(:,:,i)));
        imagesc(fullimg);axis image
        M(i) = getframe;
        writeVideo(vidObj,M(i));  
    end
    close(vidObj);
    hold off;        
       
    % save segmentation result as binary file
    seg_final = uint8(zeros(480,480,632));
    seg_final = rect_volume_insert(seg_final, app.vol.seg, app.tumor.rect, app.tumor.startZ, app.tumor.endZ);
    seg_little = downsize3D(seg_final);
    fh = fopen(['results/' strResultsLabel '_fullsize_seg.bin'],'w');
    fwrite(fh,seg_little,'uint8');
    fclose(fh);
    clear seg_final seg_little

    % write out the tmap (output seg map where each value represents LS iteration number)
    tmap_final = uint8(zeros(480,480,632));
    tmap_final = rect_volume_insert(tmap_final, app.vol.tmap, app.tumor.rect, app.tumor.startZ, app.tumor.endZ);
    tmap_little = downsize3D(tmap_final);    
    fh = fopen(['results/' strResultsLabel '_tmap.bin'],'w');
    fwrite(fh,tmap_little,'uint8');
    fclose(fh);
    clear tmap_final tmap_little
        
    disp(['Executed batch task ' itemPath ' in ' num2str(tSeg/60) ' minutes.']);