%------------------------------------------------------------------------------------
% author: Dian LU, date: 27 Mar 2019
% project: HCP100 (86 subjs), resting state
% input images: FC (in both tasks), compared between different seeds (vPre vs. dPre)
%%
clear % clean up current workspace
rmpath(genpath('/applications/spm/spm12_6906'))
rmpath(genpath('/applications/spm/spm12_7219'))
addpath('/home/dl577/spm12')
addpath('/data/dl577/scripts/spm_scripts/')

%%
dir_base = '/lustre/scratch/wbic-beta/dl577/rest_HCP/filtered_subs';
cd(dir_base)

sub_dirs = dir; % list all sub-folders
sdlist = regexp({sub_dirs.name},'\d*','match'); % return sub-folders named with numbers; variable format: {{},{},{}}
sdlist = sdlist(~cellfun(@isempty, sdlist)); % apply "isempty" function to each cell

func      = 'MNINonLinear/Results/rfMRI_REST2_LR/func';
glm2_dir  = '/data/dl577/rest_HCP/model2_pairT';
emask_dir = {'/data/dl577/masks/GM_WM_analysisMask.nii'}; % explicit mask

%% what analyses?
model_specification = 1;
model_estimation    = 1;
contrast            = 1;

%% loop through different contrast
geffect = {'pos_FC', 1; 'neg_FC', 2};

for g = 1:size(geffect, 1)
    
    con_no   = geffect{g,2};
    con_name = geffect{g,1};
    %% loop through 2 seeds
    
        disp(['Calculating for the the contrast: ' con_name '.'])
        
        model_dir = fullfile(glm2_dir,  con_name);
        
        %% -------------------MODEL SPECIFICATION------------------------------------
        if model_specification == 1
            
            clear matlabbatch
            
            mkdir(model_dir)
            % specify output dir
            matlabbatch{1}.spm.stats.factorial_design.dir = {model_dir};
            
            %% write variables and load scans
            scans = {};
            for ss = 1:86
                sub      = char(sdlist{ss});
                FuncDir  = char(fullfile(dir_base, sub, func));
                GLM1_voi1_dir = char(fullfile(FuncDir, 'model',  'GLM1', 'rsFC-vPre'));
                GLM1_voi2_dir = char(fullfile(FuncDir, 'model',  'GLM1', 'rsFC-dPre'));
                
                % loop to write compared images
                disp('Comparing between:')
                disp(fullfile(GLM1_voi1_dir, ['con_000' num2str(con_no) '.nii,1']));
                disp(fullfile(GLM1_voi2_dir, ['con_000' num2str(con_no) '.nii,1']));
                % load scans
                matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(ss).scans = {...
                    fullfile(GLM1_voi1_dir, ['con_000' num2str(con_no) '.nii,1']);
                    fullfile(GLM1_voi2_dir, ['con_000' num2str(con_no) '.nii,1'])};
            end
        end
             
        disp('Scans loaded!')
        
        matlabbatch{1}.spm.stats.factorial_design.des.pt.gmsca = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.pt.ancova = 0;
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
            
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = ['vPre>dPre-' con_name];
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
            %
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = ['dPre>vPre-' con_name];
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
            matlabbatch{1}.spm.stats.con.delete = 1;
            
            % run batch
            disp('Begin to run...');
            spm_jobman('run',matlabbatch);
            disp('Contrasts done.');
        end
    end
