% Compiles the MEX level set
disp('Compiling the level set code...');
mex matrix.c -c                
mex levelSet3D.c -c            
mex matrix.o levelSet3D.o main.c -o levelset3DC

% Tests out MEX level set on test dataset
disp('Testing the level set code with test data...');
load('test_data.mat');
tic
[seg,phi,ls_vols,tmap] = levelset3DC(double(I),double(m),100,0.25,0.9,0.5,10);
toc

disp(['Level set volume after 100 iterations is ' num2str(ls_vols(end)) ]);
disp('Level set compilation succeeded.')
