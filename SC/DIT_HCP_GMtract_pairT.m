%-----------------------------------------------------------------------
% author: Dian LU, date: 21 Mar 2019
% project: HCP100
% input images: Tracts from the seed to GM, compared between different seeds (vPre vs. dPre)
clear
base_dir = '/lustre/scratch/wbic-beta/dl577/HCP100/DTI_analysis';
cd(base_dir)
sub_dirs = dir(base_dir);
model_dir = fullfile(base_dir,'model3_pairedT_smoothed');

%% what analyses?
smooth              = 1;
model_specification = 1;
model_estimation    = 1;
contrast            = 1;

%% ----------------------SMOOTHING-----------------------------------
if smooth == 1
scans = {};
% load data
for i = 1:100
    sub_folder = sub_dirs(i+2).name;
    % loop to write compared images
    scans{end+1,1} = fullfile(sub_folder, 'T1w/Diffusion.bedpostX/vPre_DiffSpace/GM_mean_matrix2.nii,1');
    scans{end+1,1} = fullfile(sub_folder, 'T1w/Diffusion.bedpostX/dPre_DiffSpace/GM_mean_matrix2.nii,1');
end % end of the loop

matlabbatch{1}.spm.spatial.smooth.data = scans;
matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';

% run batch
disp('Begin to run smoothing...');
spm_jobman('run',matlabbatch);
end

%% -------------------MODEL SPECIFICATION------------------------------------
if model_specification == 1
    
clear matlabbatch

mkdir(model_dir)

matlabbatch{1}.spm.stats.factorial_design.dir = {model_dir};
for i = 1:100
    sub_folder = sub_dirs(i+2).name;
    % loop to write compared images
    matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(i).scans = {fullfile(sub_folder, 'T1w/Diffusion.bedpostX/vPre_DiffSpace/sGM_mean_matrix2.nii,1');
        fullfile(sub_folder, 'T1w/Diffusion.bedpostX/dPre_DiffSpace/sGM_mean_matrix2.nii,1')};
end % end of the loop
matlabbatch{1}.spm.stats.factorial_design.des.pt.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.pt.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/lustre/scratch/wbic-beta/dl577/HCP100/DTI_analysis/ROIs/GM_0.3_1.25mm_bin.nii'};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% run batch
disp('Begin to run...');
spm_jobman('run',matlabbatch);
disp('Model specified.');
end
%% -------------------MODEL ESTIMATION------------------------------------
if model_estimation == 1
clear matlabbatch
matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(model_dir,'SPM.mat')}; %Select the SPM.mat
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

% run batch
disp('Begin to run...');
spm_jobman('run',matlabbatch);
disp('Model estimated.');
end
%% ------------------CONTRAST------------------------------------
if contrast == 1
clear matlabbatch
matlabbatch{1}.spm.stats.con.spmmat = {fullfile(model_dir,'SPM.mat')};

matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'vPre>dPre';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
%
matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'dPre>vPre';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;

% run batch
disp('Begin to run...');
spm_jobman('run',matlabbatch);
disp('Contrasts done.');
end
