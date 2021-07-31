%------------------------------------------------------------------------------------
% author: Dian LU, date: 22 Mar 2019
% project: HCP100
% input images: task-state FC (in both tasks), compared between different seeds (vPre vs. dPre)
clear
base_dir = '/lustre/scratch/wbic-beta/dl577/HCP100';
cd(base_dir)
sub_dirs = dir; % list all sub-folders
sdlist = regexp({sub_dirs.name},'\d*','match'); % return sub-folders named with numbers; variable format: {{},{},{}}
sdlist = sdlist(~cellfun(@isempty, sdlist)); % apply "isempty" function to each cell
%% -------------------MODEL SPECIFICATION------------------------------------
clear matlabbatch

func1 = 'MNINonLinear/Results/tfMRI_WM_LR/fslsplit/model/GLM1_WM_v3/PPI_analyses';
func2 = 'MNINonLinear/Results/tfMRI_RELATIONAL_LR/fslsplit_scans/model/GLM1_RELATIONAL_v3/PPI_analyses';
model_dir = '/data/dl577/task_HCP/HCP100/group_results/VDpreFC_pairedT';
mkdir(model_dir)
matlabbatch{1}.spm.stats.factorial_design.dir = {model_dir};

task_code = [];
n = 0; % counter
for c = -1:2:1 % code = [-1 1] for the WM and RELATIONAL task
    for i = 1:100
        sub_folder = char(sdlist{i+2});
        
        if c == -1 % WM
            func = func1;
            voi_folders  = {'vPCC_FC', 'dPCC_FC'}; % tsFC file
        else % RELATIONAL
            func = func2;
            voi_folders  = {'aPCC_FC', 'pPCC_FC'};
        end
        %% write variables and load scans
        % if that analysis exists for this subject
        if exist([base_dir '/' sub_folder '/' func]) == 7
            disp(['Checking subj: ' sub_folder])
            n = n+1; % counter
            % loop to write compared images
            disp('Comparing between:')
            disp(fullfile(base_dir, sub_folder, func, voi_folders{1}, 'con_0003.nii,1')); % the third contrast indicate the tsFC (PPI model without interaction)
            disp(fullfile(base_dir, sub_folder, func, voi_folders{2}, 'con_0003.nii,1'));
            
            matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(n).scans = {...
                fullfile(base_dir, sub_folder, func, voi_folders{1}, 'con_0003.nii,1');
                fullfile(base_dir, sub_folder, func, voi_folders{2}, 'con_0003.nii,1')};
            % Task identity for later being regressed out
            task_code = [task_code c c];
        else
            continue
        end
    end
end
disp('Scans loaded, variables written!')

matlabbatch{1}.spm.stats.factorial_design.des.pt.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.pt.ancova = 0;
% add a covariate
matlabbatch{1}.spm.stats.factorial_design.cov.c = task_code;
matlabbatch{1}.spm.stats.factorial_design.cov.cname = 'task_identity';
matlabbatch{1}.spm.stats.factorial_design.cov.iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov.iCC = 1;
%
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
disp('Model specified.');

%% -------------------MODEL ESTIMATION------------------------------------
clear matlabbatch
matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(model_dir,'SPM.mat')}; %Select the SPM.mat
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

% run batch
disp('Begin to run...');
spm_jobman('run',matlabbatch);
disp('Model estimated.');

%% ------------------CONTRAST------------------------------------
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
