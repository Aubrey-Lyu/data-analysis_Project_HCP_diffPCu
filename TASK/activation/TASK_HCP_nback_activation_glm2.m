clear
%% Group analysis for activation studies (N-back)
%Experi = 'HCP_wm_activation';
subs1 = [100307
    100408
    101107
    101309
    101915
    103111
    103414
    103818
    105014
    105115
    106016
    108828
    110411
    111312
    111716
    113619
    113922
    114419
    115320
    116524
    117122
    118528
    118730
    118932
    120111
    122317
    122620
    123117
    123925
    124422
    125525
    126325
    127630
    127933
    128127
    128632
    129028
    130013
    130316
    131722
    133019
    133928
    135225
    135932
    136833
    138534
    139637
    140925
    144832
    146432
    147737
    148335
    148840
    149337
    149539
    149741
    151223
    151526
    151627
    153025
    154734
    156637
    159340
    160123
    161731
    162733
    163129
    176542
    178950
    188347
    189450
    190031
    192540
    196750
    198451
    199655
    201111
    211417
    211720
    212318
    214423
    221319
    239944
    245333
    280739
    298051
    366446
    414229
    499566
    654754
    672756
    751348
    792564
    856766
    857263
    899885]; % bash find function searching for subs that have con_0002.nii

subs2 = [101107
    103818
    105115
    108828
    113619
    114419
    115320
    116524
    118730
    122620
    123117
    123925
    124422
    125525
    126325
    127630
    128127
    128632
    130316
    133019
    133928
    136833
    140925
    144832
    147737
    148335
    149539
    149741
    151223
    151526
    151627
    160123
    162733
    163129
    176542
    188347
    196750
    198451
    211417
    211720
    214423
    221319
    245333
    280739
    298051
    366446
    672756
    751348
    856766
    857263];% bash find function searching for subs that have con_0008.nii
    
%% What analysis?
glm2                  = 1;
estimation2           = 1;
contrast2             = 1;
result                = 0;

%% Define subject parameters and directories
% Root Directory
dir_base    = '/lustre/scratch/wbic-beta/dl577/HCP100/';
group_dir   = '/data/dl577/task_HCP/HCP100/group_results/model';
func        = '/MNINonLinear/Results/tfMRI_WM_LR/fslsplit'; % functional scan base
model       = 'model';
glm1_task   = 'GLM1_WM_v3';
glm2_task   = 'GLM2_WM_v3';
em          = '/data/dl577/masks/GM_WM_analysisMask.nii'; % explicit mask for 2-nd level analysis

%% Initiating SPM
%spm fmri
spm_jobman('initcfg');

%% GLM2 Specification
for situation = 1:2 %situation1: only two conditions to comparel; situation2: 12 conditions to compare
    Geffect   = {'bk2-bk0_correct',     1;...
        'bk0-bk2_correct',     2;...
        'bk2-bk0_wrong',       3;...
        'bk0-bk2_wrong',       4;...
        'bk2-bk0',      5;...
        'bk0-bk2',      6;...
        'correct-wrong',           7;...
        'wrong-correct',           8;...
        'correct-wrong_bk2',     9;...
        'wrong-correct_bk2',     10;...
        'correct-wrong_bk0',     11;...
        'wrong-correct_bk0',     12};
    
    clear matlabbatch
    if situation == 1
        Geffect = Geffect(1:2,:);
        subs    = subs1;
    else
        Geffect = Geffect(3:12,:);
        subs = subs2;
    end
    
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
