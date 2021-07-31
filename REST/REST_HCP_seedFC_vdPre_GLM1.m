%======================== rest_HCP 1st-level seed-based FC =====================
% for 2 seeds: vPre, dPre, HCP100 projects
% author: Dian Lu (dl577@cam.ac.uk)
% date: 25-Mar-2019
clear % clean up current workspace
rmpath(genpath('/applications/spm/spm12_6906'))
rmpath(genpath('/applications/spm/spm12_7219'))
addpath('/home/dl577/spm12')
addpath('/data/dl577/scripts/spm_scripts/')
%% Define processing steps
glm0       = 0;
estimation = 0;
contrast   = 0;
VOI_ts     = 0;
glm1       = 1;
glm1_est   = 1;
glm1_contr = 1;

%% Root Directory
dir_base    = '/lustre/scratch/wbic-beta/dl577/rest_HCP/filtered_subs';
func        = 'MNINonLinear/Results/rfMRI_REST2_LR/func'; % functional scan base
VOI_dir     = '/data/dl577/task_HCP/HCP100/ROIs/common'; % functional rois: vPre.nii, dPre.nii
cd(dir_base)
sub_dirs = dir; % list all sub-folders
sdlist = regexp({sub_dirs.name},'\d*','match'); % return sub-folders named with numbers; variable format: {{},{},{}}
sdlist = sdlist(~cellfun(@isempty, sdlist)); % apply "isempty" function to each cell

loopin = 1%:length(sdlist);
%% Initiating SPM
%spm fmri
%spm('defaults','fmri');
spm_jobman('initcfg');

%% Define subject names
sublist = dir(fullfile(dir_base));
%--------------loop base---------------------
% Model Specification: GLM0
if glm0==1    % Model Specification: GLM0
    clear matlabbatch
    for ss = loopin
        sub      = char(sdlist{ss});
        FuncDir  = char(fullfile(dir_base, sub, func));
        GLM0_dir = char(fullfile(FuncDir, 'model', 'GLM0'));
        mkdir(GLM0_dir)
        % select functional scans once for all
        f     = spm_select('List', FuncDir, '^vol.*\.nii$');
        files = cellstr([repmat([FuncDir '/'],size(f,1),1) f repmat(',1',size(f,1),1)]);
        
        % Specify the output directory where the SPM.mat file goes.
        matlabbatch{1}.spm.stats.fmri_spec.dir = {GLM0_dir};
        % Scanner info.
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 0.72; % remember to specify
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        % load the scans
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = files;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name',{}, 'val',{});
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = cellstr(fullfile(dir_base,sub,'MNINonLinear/Results/rfMRI_REST2_LR/rfMRI_REST2_CSF_WM.txt'));
        matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
        % Specify experiment parameters
        matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{1}.spm.stats.fmri_spec.mthresh = -Inf;
        matlabbatch{1}.spm.stats.fmri_spec.mask = {'/data/dl577/masks/GM_WM_analysisMask.nii'};
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
        %
        disp('Begin to run...')
        spm_jobman('run',matlabbatch);
        disp('GLM0 specification completed.');
    end
end

%=================== Estimation, CSF/WM extraction ==================
if estimation ==1
    for ss = loopin
        sub      = char(sdlist{ss});
        FuncDir  = char(fullfile(dir_base, sub, func));
        GLM0_dir = char(fullfile(FuncDir, 'model', 'GLM0'));
        clear matlabbatch
        % Model estimation
        disp(['Estimate GLM0 for Subject: ', sub]);
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {[GLM0_dir,'/SPM.mat']}; %Select the SPM.mat
        matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        % run batch
        disp('Begin to run...');
        spm_jobman('run',matlabbatch);
        disp('GLM0 estimation job completed.');
    end
end
%% Contrast
if contrast == 1
    for ss = loopin
        sub      = char(sdlist{ss});
        FuncDir  = char(fullfile(dir_base,  sub, func));
        GLM0_dir = char(fullfile(FuncDir, 'model', 'GLM0'));       
        clear matlabbatch
        disp(['RUNNING UP Contrast for subject ' sub]);
        % define directories in the session
        matlabbatch{1}.spm.stats.con.spmmat = {[GLM0_dir,'/SPM.mat']};
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'CSF';
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'WM';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = 1;
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.delete = 1;
        % run batch
        disp('Begin to run...');
        spm_jobman('run',matlabbatch);
        disp('Contrasts job completed.');
    end
end

%% Extract timeseries of VOIs
if VOI_ts == 1
    for ss = loopin
        sub      = char(sdlist{ss});
        FuncDir  = char(fullfile(dir_base, sub, func));
        GLM0_dir = char(fullfile(FuncDir, 'model', 'GLM0'));
       
        for voi = {'vPre', 'dPre'}
            voi = char(voi);
            clear matlabbatch 
            disp(['Extract ' voi ' Timeseries  for subject ', sub '.']);
            matlabbatch{1}.spm.util.voi.spmmat = {[GLM0_dir '/SPM.mat']};
            matlabbatch{1}.spm.util.voi.adjust = NaN;
            matlabbatch{1}.spm.util.voi.session = 1;
            matlabbatch{1}.spm.util.voi.name = voi;
            matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {fullfile(VOI_dir, [voi '.nii,1'])};
            matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.5;
            matlabbatch{1}.spm.util.voi.expression = 'i1';
            % save batch file for review
            spm_jobman('run',matlabbatch);
            disp([voi ' extraction job completed.']);
        end
    end
end

%% GLM1 for seed-based FC
if glm1 == 1
    for voi = {'vPre', 'dPre'}
        voi = char(voi);
        for ss = loopin
            sub      = char(sdlist{ss});
            FuncDir  = char(fullfile(dir_base, sub, func));
            GLM1_dir = char(fullfile(FuncDir, 'model',  'GLM1', [voi 'rsFC-' voi]));
            mkdir(GLM1_dir)
            % select functional scans once for all
            f     = spm_select('List', FuncDir, '^vol.*\.nii$');
            files = cellstr([repmat([FuncDir '/'],size(f,1),1) f repmat(',1',size(f,1),1)]);
            
            clear matlabbatch
            % Specify the output directory where the SPM.mat file goes.
            matlabbatch{1}.spm.stats.fmri_spec.dir = {GLM1_dir};
            % Scanner info.
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 0.72; % remember to specify
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
            % load the scans
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = files;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress.name = voi;
            % load VOI timeseries
            load(fullfile(FuncDir, 'model', 'GLM0', ['VOI_' voi '_1.mat']))
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress.val = Y;
            clear Y
            % confounds
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = cellstr(fullfile(dir_base,sub,'MNINonLinear/Results/rfMRI_REST2_LR/rfMRI_REST2_CSF_WM.txt'));
            matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
            % Specify experiment parameters
            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = -Inf;
            matlabbatch{1}.spm.stats.fmri_spec.mask = {'/data/dl577/masks/GM_WM_analysisMask.nii'};
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
            %
            disp('Begin to run...')
            spm_jobman('run',matlabbatch);
            disp(['GLM1 specification completed for sub: ' sub ' , voi ' voi]);
        end
    end
end


%% GLM1 estimation
if glm1_est == 1
    for voi = {'vPre', 'dPre'}
        voi = char(voi);
        
        for ss = loopin
            sub      = char(sdlist{ss});
            FuncDir  = char(fullfile(dir_base, sub, func));
            GLM1_dir = char(fullfile(FuncDir, 'model',  'GLM1', [voi 'rsFC-' voi]));
            clear matlabbatch
            % Model estimation
            disp(['Estimate GLM0 for Subject: ', sub ' , voi ' voi]);
            matlabbatch{1}.spm.stats.fmri_est.spmmat = {[GLM1_dir,'/SPM.mat']}; %Select the SPM.mat
            matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
            % run batch
            disp('Begin to run...');
            spm_jobman('run',matlabbatch);
            disp('GLM1 estimation job completed.');
        end
    end
end
%% GLM1 contrast
if glm1_contr == 1
    for voi = {'vPre', 'dPre'}
        voi = char(voi);
        for ss = loopin
            sub      = char(sdlist{ss});
            FuncDir  = char(fullfile(dir_base, sub, func));
            GLM1_dir = char(fullfile(FuncDir, 'model',  'GLM1', [voi 'rsFC-' voi]));
            clear matlabbatch
            disp(['GLM1 Contrast for subject ' sub ' , voi ' voi]);
            % define directories in the session
            matlabbatch{1}.spm.stats.con.spmmat = {[GLM1_dir,'/SPM.mat']};
            %
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = voi;
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = [voi '_negative'];
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = -1;
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
            matlabbatch{1}.spm.stats.con.delete = 1;
            % run batch
            disp('Begin to run...');
            spm_jobman('run',matlabbatch);
            disp('Contrasts job completed.');
        end
    end
end
