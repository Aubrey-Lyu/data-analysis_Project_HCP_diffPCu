%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% task-HCP N-back: GLM1 activation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% version 3: modeled as event-related
clear
Experi = 'HCP_2back_';
rmpath(genpath('/applications/spm/spm12_6906'))
addpath('/home/dl577/spm12')
%% WHAT ANALYSES?
glm1     = 1;
estimate = 1;
contrast = 1;

%% Define subject parameters and directories
% Root Directory
dir_base    = '/...';
func        = '...'; % functional scan base
struc       = '...';  % structural scan base
task_model  = 'GLM1_WM_v3';
loopin      = [1:39,41:78,80:89,91:100];


%% Define subject names
sublist = dir(dir_base);
addpath /data/dl577/scripts/spm_scripts

%% Initiating SPM
%spm fmri
spm_jobman('initcfg');
%% Loop through subjects
if glm1 == 1
    for ss = loopin
        clear f files back2_corr back0_corr back2_wrong back0_wrong RT_2bk_wrong RT_2bk_corr RT_0bk_wrong RT_0bk_corr TAB block1_TRs block2_TRs block3_TRs
        % specify the current folder
        sub      = sublist(ss+2).name;
        FuncDir  = fullfile(dir_base, sub, func);
        StrucDir = fullfile(struc, sub);
        GLM1_dir = fullfile(FuncDir, 'model', task_model); %task onset with only correct trials
        tabfile = fullfile(dir_base, sub,'MNINonLinear/Results/tfMRI_WM_LR/WM_run2_TAB.txt');
        
        disp(['Specify Model for Subject: ', sub]);
        mkdir(GLM1_dir);
        
        clear matlabbatch
        if exist(tabfile,'file')==2
            % Get ready the functional scans
            cd(FuncDir); % CURRENT working directory (individual functional fmri folder)
            files = {};
            % select functional scans once for all
            f     = spm_select('List', FuncDir, '^s6vol.*\.nii$');
            files = cellstr([repmat([FuncDir '/'],size(f,1),1) f repmat(',1',size(f,1),1)]);
            
            % select experimental condition
            TAB = importfile(fullfile(dir_base, sub,'MNINonLinear/Results/tfMRI_WM_LR/WM_run2_TAB.txt'));
            back2_corr = (TAB.StimOnsetTime(~isnan(TAB.StimOnsetTime) & TAB.BlockType==2 & TAB.StimACC==1 & TAB.StimRT > 50)-TAB.SyncSlideOnsetTime(1))/1000;
            back2_wrong = (TAB.StimOnsetTime(~isnan(TAB.StimOnsetTime) & TAB.BlockType==2 & TAB.StimACC==0  & TAB.StimRT > 50)-TAB.SyncSlideOnsetTime(1))/1000;
            
            back0_corr = (TAB.StimOnsetTime(~isnan(TAB.StimOnsetTime) & TAB.BlockType==0 & TAB.StimACC==1 & TAB.StimRT > 50)-TAB.SyncSlideOnsetTime(1))/1000;
            back0_wrong = (TAB.StimOnsetTime(~isnan(TAB.StimOnsetTime) & TAB.BlockType==0 & TAB.StimACC==0 & TAB.StimRT > 50)-TAB.SyncSlideOnsetTime(1))/1000;
            
            % select session break point
            fix_onset = (TAB.Fix15secOnsetTime(strcmp(TAB.ProcedureBlock,'Fix15secPROC'))-TAB.SyncSlideOnsetTime(1))/1000;
            fix_ind = find(strcmp(TAB.ProcedureBlock,'Fix15secPROC'));
            cue_ind = fix_ind + 1;
            events = [TAB.Cue2BackOnsetTime; TAB.CueTargetOnsetTime; TAB.Fix15secOnsetTime; TAB.StimOnsetTime];
            events = ([NaN(4,1); sort(events(~isnan(events)))]-TAB.SyncSlideOnsetTime(1))/1000; % the first 4 NaN correspond to the dummy scans, so to sync with the TAB.ProcedureBlock column
            dur = (events(cue_ind(1:end-1)) - events(fix_ind(1:end-1)))/2;
            session_breaks = events(fix_ind(1:end-1)) + dur;
            tr_breaks     = round(session_breaks/0.72, 0);
            trs            = zeros(length(files),1);
            block1_TRs     = trs;
            block1_TRs(5:tr_breaks(1)) = 1;
            block2_TRs     = trs;
            block2_TRs(tr_breaks(1):tr_breaks(2)) = 1;
            block3_TRs     = trs;
            block3_TRs(tr_breaks(2):tr_breaks(3)) = 1;
            
            % Specify the output directory where the SPM.mat file goes.
            
            %% Model specification
            matlabbatch{1}.spm.stats.fmri_spec.dir = {GLM1_dir};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 0.72; % Specify the repetition time
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
            % select scans and assign to job
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = files;
            
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name ='2-back_corr';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset = back2_corr;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).orth = 1;
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
            
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name ='0-back_corr';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset = back0_corr;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).orth = 1;
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
            
            if  isempty(back2_wrong)==0 && isempty(back0_wrong)==0
                
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).name ='2-back_wrong';
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).onset = back2_wrong;
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).duration = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).tmod = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).orth = 1;
                matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
                
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).name ='0-back_wrong';
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).onset = back0_wrong;
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).duration = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).tmod = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).orth = 1;
                matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
            end
            % define regressors: WM, CSF 6 movement regressors.
            clear Y WM CSF
            % 1. specify WM.mat files
            load(fullfile(FuncDir, 'model', 'GLM0', 'VOI_WM_1.mat'))
            WM = Y;
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).name = 'WM';
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).val = WM;
            % 2. specify CSF.mat files
            load(fullfile(FuncDir, 'model', 'GLM0', 'VOI_CSF_1.mat'))
            CSF = Y;
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).name = 'CSF';
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val = CSF;
            % 3. Block/session effects
            % 3.1
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).name = 'BLOCK1';
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).val = block1_TRs;
            % 3.2
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(4).name = 'BLOCK2';
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(4).val = block2_TRs;
            % 3.3
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(5).name = 'BLOCK3';
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(5).val = block3_TRs;
            
            % 4. specify movement regressors
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {fullfile(dir_base,sub,'MNINonLinear/Results/tfMRI_WM_LR/Movement_Regressors_dt.txt')};
            % other parameters (default)
            matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = -Inf;
            matlabbatch{1}.spm.stats.fmri_spec.mask = {'/data/dl577/masks/GM_WM_analysisMask.nii'};
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
            
            % run batch
            disp('Begin to run...');
            spm_jobman('run',matlabbatch);
            disp([Experi 'glm1 specification job completed.']);
        end
    end
end

%% Model estimation
if estimate == 1
    for ss = loopin
        clear matlabbatch
        % specify the current folder
        sub      = sublist(ss+2).name;
        FuncDir  = fullfile(dir_base, sub, func);
        GLM1_dir = fullfile(FuncDir, 'model', task_model);
        tabfile = fullfile(dir_base, sub,'MNINonLinear/Results/tfMRI_WM_LR/WM_run2_TAB.txt');
        
        if exist(tabfile,'file')==2
        disp(['Estimate glm1 for Subject: ', sub]);
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {[GLM1_dir,'/SPM.mat']}; %Select the SPM.mat of glm1
        matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        
        % run batch
        disp('Begin to run...');
        spm_jobman('run',matlabbatch);
        disp([Experi 'glm1 estimation job completed.']);
        end
    end
end

%% T Contrast for First-level
if contrast == 1
    for ss = loopin
        % specify the current folder
        sub      = sublist(ss+2).name;
        FuncDir  = fullfile(dir_base, sub, func);
        GLM1_dir = fullfile(FuncDir, 'model', task_model);
        tabfile = fullfile(dir_base, sub,'MNINonLinear/Results/tfMRI_WM_LR/WM_run2_TAB.txt');
        if exist(tabfile,'file')==2
            
        clear matlabbatch
       
        matlabbatch{1}.spm.stats.con.spmmat = {[GLM1_dir,'/SPM.mat']};
        
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = '2bk-0bk(correct)';
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = '0bk-2bk(correct)';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        
        if  exist(fullfile(GLM1_dir,'beta_0022.nii'),'file')==2
                
                matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = '2bk-0bk(wrong)';
                matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [0 0 1 -1];
                matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
                
                matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = '0bk-2bk(wrong)';
                matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [0 0 -1 1];
                matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
                
                matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = '2BK-0BK';
                matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [1 -1 1 -1];
                matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
                
                matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = '0BK-2BK';
                matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = [-1 1 -1 1];
                matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
                
                matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = 'correct-wrong';
                matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = [1 1 -1 -1];
                matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
                
                matlabbatch{1}.spm.stats.con.consess{8}.tcon.name = 'wrong-correct';
                matlabbatch{1}.spm.stats.con.consess{8}.tcon.weights = [-1 -1 1 1];
                matlabbatch{1}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
                
                matlabbatch{1}.spm.stats.con.consess{9}.tcon.name = 'correct-wrong(2BK)';
                matlabbatch{1}.spm.stats.con.consess{9}.tcon.weights = [1 0 -1 0];
                matlabbatch{1}.spm.stats.con.consess{9}.tcon.sessrep = 'none';
                
                matlabbatch{1}.spm.stats.con.consess{10}.tcon.name = 'wrong-correct(2BK)';
                matlabbatch{1}.spm.stats.con.consess{10}.tcon.weights = [-1 0 1 0];
                matlabbatch{1}.spm.stats.con.consess{10}.tcon.sessrep = 'none';
                
                matlabbatch{1}.spm.stats.con.consess{11}.tcon.name = 'correct-wrong(0BK)';
                matlabbatch{1}.spm.stats.con.consess{11}.tcon.weights = [0 1 0 -1];
                matlabbatch{1}.spm.stats.con.consess{11}.tcon.sessrep = 'none';
                
                matlabbatch{1}.spm.stats.con.consess{12}.tcon.name = 'wrong-correct(0BK)';
                matlabbatch{1}.spm.stats.con.consess{12}.tcon.weights = [0 -1 0 1];
                matlabbatch{1}.spm.stats.con.consess{12}.tcon.sessrep = 'none';
                
        end
            
        matlabbatch{1}.spm.stats.con.delete = 1;
        
        % run batch
        disp('Begin to run...');
        disp([Experi ': glm1 contrast job completed.']);
        spm_jobman('run',matlabbatch);
        end
    end
end



