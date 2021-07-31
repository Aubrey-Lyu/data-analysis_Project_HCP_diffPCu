%-----------------------------------------------------------------------
% author: Dian LU, date: 21 Mar 2019
% project: HCP100
% input images: SC Tracts to GM mask(from the vPre or dPre); group analysis - 1 sample T test
clear
base_dir = '/lustre/scratch/wbic-beta/dl577/HCP100/DTI_analysis';
cd(base_dir)
sub_dirs = dir(base_dir);
glm2_dir = fullfile('/data/dl577/HCP100_DTI/model4_meanT_smoothed');
emask_dir = {'/lustre/scratch/wbic-beta/dl577/HCP100/DTI_analysis/ROIs/GM_0.3_1.25mm_bin.nii,1'}; % explicit mask
%% what analyses?
model_specification = 1;
model_estimation    = 1;
contrast            = 1;

%% loop thourgh 2 seeds
for voi = {'vPre', 'dPre'}
    voi = char(voi);
    disp(['Calculating for the seed: ' voi '.'])
    %% -------------------MODEL SPECIFICATION------------------------------------
    if model_specification == 1
        
        clear matlabbatch
        model_dir = fullfile(glm2_dir, voi);
        mkdir(model_dir)
        
        scans = {};
        for i = 1:100
            sub_folder = sub_dirs(i+2).name;
            % loop to write compared images
            scans{end+1,1} = fullfile(sub_folder, ['T1w/Diffusion.bedpostX/' voi '_DiffSpace/sGM_mean_matrix2.nii,1']);
        end % end of the loop
        
        matlabbatch{1}.spm.stats.factorial_design.dir = {model_dir};
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = scans;
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = emask_dir;
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
        %
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = voi;
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        %
        matlabbatch{1}.spm.stats.con.delete = 1;
        % run batch
        disp('Begin to run...');
        spm_jobman('run',matlabbatch);
        disp('Contrasts done.');
    end
end
