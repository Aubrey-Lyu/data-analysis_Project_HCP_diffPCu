%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% task-HCP: GLM1 relational with SESSION effect corrected %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
Experi = 'HCP_relational_task-';
%% WHAT ANALYSES?
glm1     = 1;
estimate = 1;
contrast = 1;

%% Define subject parameters and directories
% Root Directory
dir_base    = '..';
func        = '...'; % functional scan base
struc       = '..';  % structural scan base
task_model  = 'GLM1_RELATIONAL_v3';
loopin      = 1:100;

%% Define subject names
sublist = dir(dir_base);
addpath /data/dl577/scripts/spm_scripts

%% Initiating SPM
%spm fmri
spm_jobman('initcfg');
%% Loop through subjects
if glm1 == 1
    for ss = loopin
        
        % specify the current folder
        sub      = sublist(ss+2).name;
        FuncDir  = fullfile(dir_base, sub,  func);
        StrucDir = fullfile(struc, sub);
        GLM1_dir = fullfile(FuncDir, 'model', task_model); %task onset with only correct trials
        tabfile = fullfile(dir_base, sub,'MNINonLinear/Results/tfMRI_RELATIONAL_LR//RELATIONAL_cleaned_TAB2.txt'); % with behavioural response information
        clear matlabbatch f TAB
        
        if exist(tabfile,'file')==2
            
            % Get ready the functional scans
            cd(FuncDir); % CURRENT working directory (individual functional fmri folder)
            files = {};
            % select functional scans once for all
            f     = spm_select('List', FuncDir, '^s6vol0.*\.nii$');
            files = cellstr([repmat([FuncDir '/'],size(f,1),1) f repmat(',1',size(f,1),1)]);
            
            % select experimental condition
            TAB = importfile_relat(tabfile, 2,37);
            relat1 = (TAB.rela_onset(strcmp(TAB.Procedure, 'RelationalPROC') & TAB.rela_acc==1& TAB.RelationalSlideRT~=0))/1000;
            relat0 = (TAB.rela_onset(strcmp(TAB.Procedure,'RelationalPROC') & TAB.rela_acc==-1& TAB.RelationalSlideRT~=0))/1000;
            contr1 = (TAB.cont_onset(strcmp(TAB.Procedure,'ControlPROC') & TAB.cont_acc==1& TAB.RelationalSlideRT~=0))/1000;
            contr0 = (TAB.cont_onset(strcmp(TAB.Procedure,'ControlPROC') & TAB.cont_acc==-1& TAB.RelationalSlideRT~=0))/1000;
            dur_relat1 = (TAB.RelationalSlideRT(strcmp(TAB.Procedure,'RelationalPROC') & TAB.rela_acc==1 & TAB.RelationalSlideRT~=0))/1000;
            dur_relat0 = (TAB.RelationalSlideRT(strcmp(TAB.Procedure,'RelationalPROC') & TAB.rela_acc==-1& TAB.RelationalSlideRT~=0))/1000;
            dur_contr1 = (TAB.ControlSlideRT(strcmp(TAB.Procedure,'ControlPROC') & TAB.cont_acc==1& TAB.RelationalSlideRT~=0))/1000;
            dur_contr0 = (TAB.ControlSlideRT(strcmp(TAB.Procedure,'ControlPROC') & TAB.cont_acc==-1& TAB.RelationalSlideRT~=0))/1000;
            dur_fix    = nanmean(TAB.fix_offset-TAB.fix_onset)/2;
            session_breaks = TAB.fix_onset + dur_fix;
            session_breaks = session_breaks(~isnan(session_breaks));
            tr_breaks      = round(session_breaks/720, 0);  

            trs            = zeros(length(files),1);
            block1_TRs     = trs;
            block1_TRs(5:tr_breaks(1)) = 1;
            block2_TRs     = trs;
            block2_TRs(tr_breaks(1):tr_breaks(2)) = 1;
            block3_TRs     = trs;
            block3_TRs(tr_breaks(2):tr_breaks(3)) = 1;
            
            % Specify the output directory where the SPM.mat file goes.
            
            %% Model specification
                        
            disp(['Specify Model for Subject: ', sub]);
            mkdir(GLM1_dir);
            
            matlabbatch{1}.spm.stats.fmri_spec.dir = {GLM1_dir};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 0.72; % Specify the repetition time
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
            % select scans and assign to job
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = files;
            % select experimental condition
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name ='relational-correct';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset = relat1;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = dur_relat1;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).orth = 1;
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
            
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name ='control-correct';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset = contr1;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = dur_contr1;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).orth = 1;
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
            if  isempty(relat0)==0 && isempty(contr0)==0
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).name ='relational-wrong';
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).onset = relat0;
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).duration = dur_relat0;
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).tmod = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).orth = 1;
                matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
                
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).name ='control-wrong';
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).onset = contr0;
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).duration = dur_contr0;
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).tmod = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).orth = 1;
                matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
            end
            % define regressors: WM, CSF 6 movement regressors.
            %
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
            %
            % 4. specify movement regressors
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {fullfile(dir_base,sub,'MNINonLinear/Results/tfMRI_RELATIONAL_LR/Movement_Regressors.txt')};
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
        FuncDir  = fullfile(dir_base, sub,func);
        GLM1_dir = fullfile(FuncDir, 'model', task_model);
        tabfile = fullfile(dir_base, sub,'MNINonLinear/Results/tfMRI_RELATIONAL_LR/RELATIONAL_cleaned_TAB2.txt');
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
        FuncDir  = fullfile(dir_base, sub,  func);
        GLM1_dir = fullfile(FuncDir, 'model', task_model);
        tabfile = fullfile(dir_base, sub,'MNINonLinear/Results/tfMRI_RELATIONAL_LR/RELATIONAL_cleaned_TAB2.txt');
        if exist(tabfile,'file')==2
            clear matlabbatch
            matlabbatch{1}.spm.stats.con.spmmat = {[GLM1_dir,'/SPM.mat']};
            
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'relational-control(correct)';
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
            
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'control-relational(correct)';
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
            
            if  exist(fullfile(GLM1_dir,'beta_0022.nii'),'file')==2
                
                matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'relational-control(wrong)';
                matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [0 0 1 -1];
                matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
                
                matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'control-relational(wrong)';
                matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [0 0 -1 1];
                matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
                
                matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'relational-control';
                matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [1 -1 1 -1];
                matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
                
                matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'control-relational';
                matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = [-1 1 -1 1];
                matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
                
                matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = 'correct-wrong';
                matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = [1 1 -1 -1];
                matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
                
                matlabbatch{1}.spm.stats.con.consess{8}.tcon.name = 'wrong-correct';
                matlabbatch{1}.spm.stats.con.consess{8}.tcon.weights = [-1 -1 1 1];
                matlabbatch{1}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
                
                matlabbatch{1}.spm.stats.con.consess{9}.tcon.name = 'correct-wrong(relational)';
                matlabbatch{1}.spm.stats.con.consess{9}.tcon.weights = [1 0 -1 0];
                matlabbatch{1}.spm.stats.con.consess{9}.tcon.sessrep = 'none';
                
                matlabbatch{1}.spm.stats.con.consess{10}.tcon.name = 'wrong-correct(relational)';
                matlabbatch{1}.spm.stats.con.consess{10}.tcon.weights = [-1 0 1 0];
                matlabbatch{1}.spm.stats.con.consess{10}.tcon.sessrep = 'none';
                
                matlabbatch{1}.spm.stats.con.consess{11}.tcon.name = 'correct-wrong(control)';
                matlabbatch{1}.spm.stats.con.consess{11}.tcon.weights = [0 1 0 -1];
                matlabbatch{1}.spm.stats.con.consess{11}.tcon.sessrep = 'none';
                
                matlabbatch{1}.spm.stats.con.consess{12}.tcon.name = 'wrong-correct(control)';
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



