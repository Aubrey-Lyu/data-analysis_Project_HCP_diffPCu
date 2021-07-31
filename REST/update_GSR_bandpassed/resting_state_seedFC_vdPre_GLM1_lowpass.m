%======================== rest_HCP 1st-level seed-based FC =====================
% for 2 seeds: vPre, dPre, HCP100 projects
% author: Dian Lu (dl577@cam.ac.uk)
% date: 25-Mar-2019
clear % clean up current workspace
tic
addpath('/home/dl577/spm12')
addpath(genpath('/home/dl577/scripts/spm_scripts/'))
%% Define processing steps
glm0       = 1;
estimation = 1;
contrast   = 1;
VOI_ts     = 1;
glm1       = 1;
glm1_est   = 1;
glm1_contr = 1;

%% Root Directory
dir_base    = '/rds/project/rds-pXaBn8E6hyM/migration/group/uda/ccig/dl577/rest_HCP/filtered_subs';
func        = 'MNINonLinear/Results/rfMRI_REST2_LR/func'; % functional scan base
VOI_dir     = '/rds/user/dl577/rds-uda-2-pXaBn8E6hyM/migration/data/dl577/task_HCP/HCP100/ROIs/common'; % functional rois: vPre.nii, dPre.nii
BrainMask_dir = '/rds/user/dl577/rds-uda-2-pXaBn8E6hyM/migration/data/dl577/whole-brain-mask';
batch_folder = '/rds/project/rds-pXaBn8E6hyM/migration/group/uda/ccig/dl577/rest_HCP/matlabbatches';

cd(dir_base)
%% Define subject names
sub_dirs = dir; % list all sub-folders
sdlist = regexp({sub_dirs.name},'\d*','match'); % return sub-folders named with numbers; variable format: {{},{},{}}
sdlist = sdlist(~cellfun(@isempty, sdlist)); % apply "isempty" function to each cell

loopin = 1:length(sdlist);

%% Initiating SPM
%spm fmri
%spm('defaults','fmri');
spm_jobman('initcfg');

%--------------loop base---------------------
%% Model Specification: GLM0
if glm0==1    % Model Specification: GLM0
    
    counter = 0;
    matlabbatch_glm0 = {};
    for ss = loopin
        sub      = char(sdlist{ss});
        FuncDir  = char(fullfile(dir_base, sub, func));
        GLM0_dir = char(fullfile(FuncDir, 'model', 'GLM0_LP'));
        mkdir(GLM0_dir)
        
        if exist(fullfile(GLM0_dir, 'SPM.mat'), 'file') ~= 2
            counter = counter +1;
          clear matlabbatch
           
            % select functional scans once for all
            files = cellstr(expand_4d_vols(fullfile(FuncDir, 'brikFile.nii'))); % this is low_pass data (4D file)
            %f     = spm_select('List', FuncDir, '^vol.*\.lowpassed.nii$');
            %files = cellstr([repmat([FuncDir '/'],size(f,1),1) f repmat(',1',size(f,1),1)]); % (3D file)
            
            %%
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
            %matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = cellstr(fullfile(dir_base,sub,'MNINonLinear/Results/rfMRI_REST2_LR/rfMRI_REST2_CSF_WM.txt'));
            matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128; % high pass filter
            
            % Specify experiment parameters
            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = -Inf;
            matlabbatch{1}.spm.stats.fmri_spec.mask = {fullfile(BrainMask_dir, 'GM_WM_mask.nii')};
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
             
                % -------store subj-specific parameters to a meta-batch for parallelising job ----------
            
            matlabbatch_glm0{counter} = matlabbatch;
            clear matlabbatch
            
        end
    end
     % save batch file for review
    save(fullfile(batch_folder, 'matlabbatch_glm0.mat'), 'matlabbatch_glm0')
    
    % ===== run file ========
    batch = matlabbatch_glm0;
    parfor counter = 1: length(matlabbatch_glm0) % parallelising the job
        disp(['Running for GLM0: batch-' num2str(counter)])
        spm_jobman('run', batch{counter});
    end
    disp('GLM0 specification completed.');
    clear batch
end

%% ==================== Estimation ===================
if estimation ==1
    counter = 0;
    matlabbatch_est0 = {};
    clear matlabbatch
    for ss = loopin
        sub      = char(sdlist{ss});
        FuncDir  = char(fullfile(dir_base, sub, func));
        GLM0_dir = char(fullfile(FuncDir, 'model', 'GLM0_LP'));
        
        if exist(fullfile(GLM0_dir, 'ResMS.nii'), 'file') ~= 2
            counter = counter +1;
            % Model estimation
            disp(['Estimate GLM0 for Subject: ', sub]);
            matlabbatch{1}.spm.stats.fmri_est.spmmat = {[GLM0_dir '/SPM.mat']}; %Select the SPM.mat
            matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
            % -------store subj-specific parameters to a meta-batch for parallelising job ----------
            
            matlabbatch_est0{counter} = matlabbatch;
            clear matlabbatch
        end
    end
     % save batch file for review
    save(fullfile(batch_folder, 'matlabbatch_est0.mat'), 'matlabbatch_est0')
    
    %% ===== run batch ========
    batch = matlabbatch_est0;
    parfor counter = 1: length(batch)
        disp(['Running for GLM0 estimation: batch-' num2str(counter)])
        spm_jobman('run', batch{counter});
    end
    clear batch
    disp('GLM0 estimation job completed.');
end

%% ===== Contrast ======
if contrast == 1
    counter = 0;
    matlabbatch_ctr0 = {};
    for ss = loopin
        sub      = char(sdlist{ss});
        FuncDir  = char(fullfile(dir_base,  sub, func));
        GLM0_dir = char(fullfile(FuncDir, 'model', 'GLM0_LP'));
        %------------------------------------------------------------%
        if exist(fullfile(GLM0_dir, 'con_0002.nii'), 'file') ~= 2
            counter = counter + 1;
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
             
                % -------store subj-specific parameters to a meta-batch for parallelising job ----------
            matlabbatch_ctr0{counter} = matlabbatch;
            clear matlabbatch
        end
    end
     % save batch file for review
    save(fullfile(batch_folder, 'matlabbatch_ctr0.mat'), 'matlabbatch_ctr0')
    %% run batch
    disp('Begin to run...');
    batch = matlabbatch_ctr0;
    parfor counter = 1: length(batch)
        disp(['Running for GLM0 contrast: batch-' num2str(counter)])
        spm_jobman('run',batch{counter});
    end
    clear batch
    disp('Contrasts job completed.');
end

%% ===== Extract timeseries of VOIs, (implemented global signal extraction)=====
if VOI_ts == 1
    matlabbatch_voi = {};
    counter = 0;
    for ss = loopin
        sub      = char(sdlist{ss});
        FuncDir  = char(fullfile(dir_base, sub, func));
        GLM0_dir = char(fullfile(FuncDir, 'model', 'GLM0_LP'));
        
        for voi = {'vPre', 'dPre', 'TPM0000', 'TPM0001', 'TPM0002'}
            voi = char(voi);
            %------------------------------------------------------------%
            if exist(fullfile(GLM0_dir, ['VOI_' voi '_1.mat']), 'file') ~= 2
                counter = counter + 1;
                clear matlabbatch
                %
                disp(['Extract ' voi ' Timeseries  for subject ', sub '.']);
                matlabbatch{1}.spm.util.voi.spmmat = {[GLM0_dir '/SPM.mat']};
                matlabbatch{1}.spm.util.voi.adjust = NaN; % adjust for everything
                matlabbatch{1}.spm.util.voi.session = 1;
                matlabbatch{1}.spm.util.voi.name = voi;
                if ismember(voi, {'vPre', 'dPre'})
                    matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {fullfile(VOI_dir, [voi '.nii,1'])};
                else
                    matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {fullfile(BrainMask_dir, [voi '.nii,1'])};
                end
                matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.5;
                matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {fullfile(GLM0_dir,  'mask.nii,1')};
                matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5;
                %
                matlabbatch{1}.spm.util.voi.expression = 'i1&i2';
                
                % -------store subj-specific parameters to a meta-batch for parallelising job ----------
                 matlabbatch_voi{counter} = matlabbatch;
            end
        end
    end
     % save batch file for review
    save(fullfile(batch_folder, 'matlabbatch_voi.mat'), 'matlabbatch_voi')
    
    %% ========== run batch =================
    disp('Begin to run...');
    batch = matlabbatch_voi;
    parfor counter = 1: length(batch)
        disp(['Running for VOI ts extraction: batch-' num2str(counter)])
        spm_jobman('run',batch{counter});
    end
    clear batch
    disp('Time-series extraction job completed.');
end

%% ============== GLM1 for seed-based FC =================
if glm1 == 1
    matlabbatch_glm1={};
    counter = 0;
    for voi = {'vPre', 'dPre'}
        voi = char(voi);
        for ss = loopin
            sub      = char(sdlist{ss});
            FuncDir  = char(fullfile(dir_base, sub, func));
            GLM1_dir = char(fullfile(FuncDir, 'model',  'GLM1_LP', ['rsFC-' voi]));
            mkdir(GLM1_dir)
            
            if exist(fullfile(GLM1_dir, 'SPM.mat'), 'file') ~= 2
                % ------------------------
                counter = counter +1;
                % select functional scans once for all
                files = cellstr(expand_4d_vols(fullfile(FuncDir, 'brikFile.nii')));
                % f     = spm_select('List', FuncDir, '^vol.*\.lowpassed.nii$');
                % files = cellstr([repmat([FuncDir '/'],size(f,1),1) f repmat(',1',size(f,1),1)]);
                
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
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).name = voi;
                % load VOI timeseries
                load(fullfile(FuncDir, 'model', 'GLM0_LP', ['VOI_' voi '_1.mat']))
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).val = Y;
                clear Y
                
                % gloabl signal confound
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).name = 'GM';
                % load confound  voi timeseries
                load(fullfile(FuncDir, 'model', 'GLM0_LP', 'VOI_TPM0000_1.mat'))
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val = Y;
                clear Y
                
                % white matter
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).name = 'WM';
                % load confound  voi timeseries
                load(fullfile(FuncDir, 'model', 'GLM0_LP', 'VOI_TPM0001_1.mat'))
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).val = Y;
                clear Y
                
                % CSF
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(4).name = 'CSF';
                % load confound  voi timeseries
                load(fullfile(FuncDir, 'model', 'GLM0_LP', 'VOI_TPM0002_1.mat'))
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(4).val = Y;
                clear Y
                % confounds
                %matlabbatch{counter}.spm.stats.fmri_spec.sess.multi_reg = cellstr(fullfile(dir_base,sub,'MNINonLinear/Results/rfMRI_REST2_LR/rfMRI_REST2_CSF_WM.txt'));
                matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
                % Specify experiment parameters
                matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
                matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
                matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
                matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
                matlabbatch{1}.spm.stats.fmri_spec.mthresh = -Inf;
                matlabbatch{1}.spm.stats.fmri_spec.mask = {fullfile(BrainMask_dir, 'GM_WM_mask.nii')};
                matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
                % -------store subj-specific parameters to a meta-batch for parallelising job ----------
                matlabbatch_glm1{counter} = matlabbatch;
                clear matlabbatch
            end
        end
    end
     % save batch file for review
    save(fullfile(batch_folder, 'matlabbatch_glm1.mat'), 'matlabbatch_glm1')
    
    %% ============= run batch ==================
    disp('Begin to run...');
    batch = matlabbatch_glm1;
    parfor counter = 1: length(batch)
        disp(['Running for GLM1-seedbased FC: batch-' num2str(counter)])
        spm_jobman('run',batch{counter});
    end
    clear batch
    disp('GLM1 job completed.');
end

%% GLM1 estimation
if glm1_est == 1
    matlabbatch_est1={};
    counter = 0;
    for voi = {'vPre', 'dPre'}
        voi = char(voi);
        
        for ss = loopin
            sub      = char(sdlist{ss});
            FuncDir  = char(fullfile(dir_base, sub, func));
            GLM1_dir = char(fullfile(FuncDir, 'model',  'GLM1_LP', ['rsFC-' voi]));
            
            if exist(fullfile(GLM1_dir, 'ResMS.nii'), 'file') ~= 2
                % -----------------------
                counter = counter +1;
                % Model estimation
                disp(['Estimate GLM0 for Subject: ', sub ' , voi ' voi]);
                matlabbatch{1}.spm.stats.fmri_est.spmmat = {[GLM1_dir,'/SPM.mat']}; %Select the SPM.mat
                matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
                matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
                
                % -------store subj-specific parameters to a meta-batch for parallelising job ----------
                matlabbatch_est1{counter} = matlabbatch;
                clear matlabbatch
            end
        end
    end
     % save batch file for review
    save(fullfile(batch_folder, 'matlabbatch_est1.mat'), 'matlabbatch_est1')
    
   %% ============= run batch ==================
    disp('Begin to run...');
    batch = matlabbatch_est1;
    parfor counter = 1: length(batch)
        disp(['Running estimation for GLM1-seedbased FC: batch-' num2str(counter)])
        spm_jobman('run', batch{counter});
    end
    disp('GLM1 estimation job completed.');
end

%% GLM1 contrast
if glm1_contr == 1
    matlabbatch_ctr1 = {};
    counter = 0;
    
    for voi = {'vPre', 'dPre'}
        voi = char(voi);
        for ss = loopin
            sub      = char(sdlist{ss});
            FuncDir  = char(fullfile(dir_base, sub, func));
            GLM1_dir = char(fullfile(FuncDir, 'model',  'GLM1_LP', ['rsFC-' voi]));
            
            if exist(fullfile(GLM1_dir, 'con_0002.nii'), 'file') ~= 2
                counter = counter +1;
                %-------------------------------
                disp(['GLM1 Contrast for subject ' sub ' , voi ' voi]);
                % define directories in the session
                matlabbatch{1}.spm.stats.con.spmmat = {[GLM1_dir,'/SPM.mat']};
                %
                matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = [voi '-FC'];
                matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
                matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
                matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = [voi '-FC_negative'];
                matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = -1;
                matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
                matlabbatch{1}.spm.stats.con.delete = 1;
                 
                % -------store subj-specific parameters to a meta-batch for parallelising job ----------
                matlabbatch_ctr1{counter} = matlabbatch;
                clear matlabbatch
            end
        end
    end
     % save batch file for review
    save(fullfile(batch_folder, 'matlabbatch_ctr1.mat'), 'matlabbatch_ctr1')
    
    %% ============= run batch ==================
    disp('Begin to run...');
    batch = matlabbatch_ctr1;
    parfor counter = 1: length(batch)
        disp(['Running contrast for GLM1-seedbased FC: batch-' num2str(counter)])
        spm_jobman('run', batch{counter});
    end
    clear batch
    disp('GLM1 contrast job completed.');
end
toc