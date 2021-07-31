%-----------------------------------------------------------------------
% second level: difference in RT as covariates
%-----------------------------------------------------------------------
clear

%% what analyses?
model_speci = 1;
estimate_contrast = 1;

%% import variables
addpath /data/dl577/scripts/spm_scripts
metadcm = import_metadcm('/lustre/scratch/wbic-beta/dl577/HCP100/metadcm.txt');
df = read_rtdiff('/lustre/scratch/wbic-beta/dl577/HCP100/rt_diff.txt');

%% Define directories
dir_base  = '/lustre/scratch/wbic-beta/dl577/HCP100/';
group_base = '/data/dl577/task_HCP/HCP100/group_results/Regression_diffRT_RFX';
tasks      = {'WM','RELATIONAL'};

Geffect = {'difficult-easy',     1;...
    'easy-difficult',     2;...
    'vPre_difficult-easy', 1;...
    'vPre_easy-difficult', 2;...
    'dPre_difficult-easy', 1;...
    'dPre_easy-difficult', 2};

for task = tasks
    task = task{1};
    for gg = 1: size(Geffect, 1)
        
        effect = Geffect{gg,1};
        c      = Geffect{gg,2};
        
        confile = ['con_000' num2str(c) '.nii,1'];
        
        output_dir = fullfile(group_base, effect,task);
        
        scans    = {};
        age_codes = [];
        gender_codes = [];
        rt_diff = [];
        
        if model_speci == 1
            
            mkdir(output_dir)
            
            if isempty(strfind(Geffect{gg,1}, 'Pre'))
                disp('Activation results...')
                func1 = 'MNINonLinear/Results/tfMRI_WM_LR/fslsplit/model/GLM1_WM_v3';
                func2 = 'MNINonLinear/Results/tfMRI_RELATIONAL_LR/fslsplit_scans/model/GLM1_RELATIONAL_v3';
                
            else
                if ~isempty(strfind(Geffect{gg,1}, 'dPre'))
                    
                    disp('PPI results of dPre...')
                    func1 = 'MNINonLinear/Results/tfMRI_WM_LR/fslsplit/model/GLM1_WM_v3/PPI_analyses/dPCC_interaction';
                    func2 =  'MNINonLinear/Results/tfMRI_RELATIONAL_LR/fslsplit_scans/model/GLM1_RELATIONAL_v3/PPI_analyses/pPCC_interaction';
                    
                    
                elseif ~isempty(strfind(Geffect{gg,1}, 'vPre'))
                    disp('PPI results of vPre...')
                    func1 = 'MNINonLinear/Results/tfMRI_WM_LR/fslsplit/model/GLM1_WM_v3/PPI_analyses/vPCC_interaction';
                    func2 =  'MNINonLinear/Results/tfMRI_RELATIONAL_LR/fslsplit_scans/model/GLM1_RELATIONAL_v3/PPI_analyses/aPCC_interaction';
                    
                end
            end
            
            % prepare variables
            switch task
                        case 'WM'
                            func = func1;
                        case 'RELATIONAL'
                            func = func2;
            end         
              
            loopin = 1 : length(df.SUBJ);
            for ss = loopin
                sub = num2str(df.SUBJ(ss));
                if strcmp(df.TASK(ss), task)
                    scan = fullfile(dir_base, sub, func, confile);
                    a = df.age_code(ss);
                    g = df.gender_code(ss);
                    r = df.RT_diff(ss);
                    
                    scans{end+1,1} = scan;
                    age_codes = [age_codes, a];
                    gender_codes = [gender_codes, g];
                    rt_diff = [rt_diff, r];
                end
            end
            
            %% create matlabbatch
            clear matlabbatch
            
            matlabbatch{1}.spm.stats.factorial_design.dir = {output_dir};
            % con_000* scans
            matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = scans;
            
            % main regressor: RT diff
            matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.c = rt_diff;
            matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.cname = 'RT_diff';
            matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 1;
            matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
            % age
            matlabbatch{1}.spm.stats.factorial_design.cov(1).c = age_codes;
            matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'AGE';
            matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
            matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
            % gender
            matlabbatch{1}.spm.stats.factorial_design.cov(2).c = gender_codes;
            matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'GENDER';
            matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
            matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 1;
            % other parameters, default
            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
            
            spm_jobman('run',matlabbatch);
        end
        
        %% ---- model estimation -------
        
        if estimate_contrast  == 1
            clear matlabbatch
            % model estimation
            matlabbatch{1}.spm.stats.fmri_est.spmmat = {[output_dir,'/SPM.mat']}; %Select the SPM.mat of glm1
            matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
            
            disp('Model contrast...')
            % model contrast
            matlabbatch{2}.spm.stats.con.spmmat = {[output_dir,'/SPM.mat']};
            
            matlabbatch{2}.spm.stats.con.consess{1}.tcon.name = 'RT_diff';
            matlabbatch{2}.spm.stats.con.consess{1}.tcon.weights = [0 0 0 1];
            matlabbatch{2}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
            
            matlabbatch{2}.spm.stats.con.consess{2}.tcon.name = 'RT_diff_nega';
            matlabbatch{2}.spm.stats.con.consess{2}.tcon.weights = [0 0 0 -1];
            matlabbatch{2}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
            
            matlabbatch{2}.spm.stats.con.delete = 0;
            
            spm_jobman('run',matlabbatch);
        end
        
        
    end
end