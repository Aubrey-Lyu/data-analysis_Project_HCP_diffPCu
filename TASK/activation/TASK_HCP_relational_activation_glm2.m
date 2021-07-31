clear
%% Group analysis for activation studies (relational processing)
%Experi = 'HCP_relational_activation';

%% What analysis?
glm2                  = 1;
estimation2           = 1;
contrast2             = 1;
result                = 0;

%% Define subject parameters and directories
% Root Directory
dir_base    = '/...';
group_dir   = '/...';
func        = '/...'; % functional scan base
model       = 'model';
glm1_task   = 'GLM1_RELATIONAL_v3';
glm2_task   = 'GLM2_RELATIONAL_v4';
em          = '../GM_WM_analysisMask.nii'; % explicit mask for 2-nd level analysis
loopin      = 1:100;
%% Define subject names
%sublist = dir(dir_base);
subs=[...];

%% Initiating SPM
%spm fmri
spm_jobman('initcfg');

%% GLM2 Specification
for situation = 1:2 %situation1: only two conditions to comparel; situation2: 12 conditions to compare
    Geffect   = {'relat-contr_correct',     1;...
                 'contr-relat_correct',     2;...
                 'relat-contr_wrong',       3;...
                 'contr-relat_wrong',       4;...
                 'relational-control',      5;...
                 'control-relational',      6;...
                 'correct-wrong',           7;...
                 'wrong-correct',           8;...
                 'correct-wrong_relat',     9;...
                 'wrong-correct_relat',     10;...
                 'correct-wrong_contr',     11;...
                 'wrong-correct_contr',     12};

    clear matlabbatch
    
    
    if glm2 == 1
        
        for i = 1:size(Geffect,1)
            gg = Geffect{i,1};
            cc = Geffect{i,2};
            if cc < 10
                cfile = ['con_000' num2str(Geffect{i,2}) '.nii,1'];
            else
                cfile = ['con_00' num2str(Geffect{i,2}) '.nii,1'];
            end
            disp(['Group-level specification for the ' gg ' effect.']);
            
            cfiles = {};
            for ss = 1:length(subs)
                sub = subs(ss);
                FuncDir  = fullfile(dir_base, num2str(sub), func);
                GLM1_dir = fullfile(FuncDir, model, glm1_task);
                cfiles{end+1,1} = fullfile(GLM1_dir, cfile);
            end
            
            GLM2_dir = fullfile(group_dir, glm2_task, 'activation_results', gg);
            mkdir(GLM2_dir);
            
            clear matlabbatch
            matlabbatch{1}.spm.stats.factorial_design.dir = {GLM2_dir};
            matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = cfiles;
            matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
            
            % run batch
            disp('Begin to run...');
            spm_jobman('run',matlabbatch);
            disp('glm2 specification completed.');
        end
    end
    clear gg i cc
    
    %% Estimation
    if estimation2 == 1

        for i = 1:size(Geffect,1)
            gg = Geffect{i,1};
            
            disp(['GLM2 estimation for the ' gg ' effect.']);
            
            clear matlabbatch
            % specify the model folder
            GLM2_dir = fullfile(group_dir, glm2_task, 'activation_results', gg);
            disp(['Estimate glm2 for the effect: ', gg]);
            matlabbatch{1}.spm.stats.fmri_est.spmmat = {[GLM2_dir,'/SPM.mat']}; %Select the SPM.mat of glm2
            matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
            
            % run batch
            disp('Begin to run...');
            spm_jobman('run',matlabbatch);
            disp('glm2 estimation job completed.');
        end
    end
    
    
    %% UP contrast to see the effect of interest
    if contrast2 == 1
        for i = 1: size(Geffect,1)
            cc = Geffect{i,2};
            gg = Geffect{i,1};
            GLM2_dir = fullfile(group_dir, glm2_task, 'activation_results', gg);
            disp(['Contrast2 for the ' gg ' effect--Contrast ' num2str(cc)]);
            clear matlabbatch
            matlabbatch{1}.spm.stats.con.spmmat = {[GLM2_dir,'/SPM.mat']};
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = gg;
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
            % run batch
            disp('Begin to run...');
            spm_jobman('run',matlabbatch);
            disp('glm2 result completed.');
            spm_jobman('run',matlabbatch);
        end
    end
    
    %% Result reports
    if result == 1
        for i = 1:size(Geffect,1)
            gg = Geffect{i,1};
            cc = Geffect{i,2};
            
            GLM2_dir = fullfile(group_dir, glm2_task, 'activation_results', gg);
            disp(['View result for the ' gg ' effect--Contrast ' num2str(cc)]);
            
            clear matlabbatch
            matlabbatch{1}.spm.stats.results.spmmat = {[GLM2_dir '/SPM.mat']};
            matlabbatch{1}.spm.stats.results.conspec.titlestr = gg;
            matlabbatch{1}.spm.stats.results.conspec.contrasts = 1;
            matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
            matlabbatch{1}.spm.stats.results.conspec.thresh = 0.001;
            matlabbatch{1}.spm.stats.results.conspec.extent = 0;
            matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
            matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
            matlabbatch{1}.spm.stats.results.units = 1;
            matlabbatch{1}.spm.stats.results.export{1}.ps = true;
            
            % run batch
            disp('Begin to run...');
            spm_jobman('run',matlabbatch);
            disp('glm2 result completed.');
            spm_jobman('run',matlabbatch);
        end
    end
    
end
