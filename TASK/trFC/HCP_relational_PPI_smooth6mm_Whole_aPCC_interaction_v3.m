clear
rmpath(genpath('/applications/spm/spm12_6906'))
addpath('/home/dl577/spm12')
Experi = 'taskRELATIONAL_PPI_6mm_aPCC_interaction_v3_';
%% What analysis?
PPI_specification     = 1;
PPI_GLM_specification = 1;
estimation            = 1;
contrast              = 1;
glm2                  = 1;
estimation2           = 1;
contrast2             = 1;
result                = 0;

%% Define subjects to loop through
sublist  = [188347
    117122
    122620
    189450
    135932
    149741
    108828
    100307
    138534
    103414
    147737
    135225
    366446
    201111
    672756
    144832
    499566
    105014
    105115
    139637
    239944
    856766
    149337
    111716
    100408
    397760
    115320
    133928
    156637
    153025
    133019
    178950
    101309
    131217
    128127
    214423
    113619
    118932
    123117
    654754
    857263
    122317
    110411
    127630
    163129
    212318
    111312
    151627
    211720
    118528
    792564
    128632
    114419
    118730
    161731
    159340
    148840
    131722
    136833
    198451
    245333
    130316
    149539
    129028
    751348
    221319
    124422
    199655
    146432
    113922
    151223
    176542
    106016
    101915
    125525
    116524
    190031
    120111
    211417
    151526
    899885
    103818
    160123
    148335
    127933
    414229
    298051
    123925
    101107
    280739
    130013
    126325
    103111
    162733
    192540
    208226
    154734
    756055];

loopin = 1:length(sublist);

%% Define directories
dir_base  = '/lustre/scratch/wbic-beta/dl577/HCP100/';
group_base= '/data/dl577/task_HCP/RELATIONAL/group_results';
model     = 'model';
func      = 'MNINonLinear/Results/tfMRI_RELATIONAL_LR/fslsplit_scans'; % functional scan base
struc     = 'T1w/Results/tfMRI_WM_LR';  % structural scan base
glm1_task = 'GLM1_RELATIONAL_v3';
glm2_task = 'GLM2_RELATIONAL_v3';
voi       = 'aPCC';
PPI_folder = 'PPI_analyses';
PPI_name  = [voi '_interaction'];
Geffect   = {'apcc_relat-contr',1;...
            'apcc_contr-relat', 2};

%% Initiating SPM
%spm fmri
%spm_jobman('initcfg');

%% Loop through subjects
if PPI_specification ==1
    for ss = loopin
        clear matlabbatch
        % specify the current folder
        sub      = num2str(sublist(ss));
        FuncDir  = fullfile(dir_base, sub, func);
        GLM1_dir = fullfile(FuncDir, model, glm1_task);
        
        disp(['PPI specification for the No.' num2str(ss) ' subject: ' sub '.'])

        matlabbatch{1}.spm.stats.ppi.spmmat = {fullfile(GLM1_dir, 'SPM.mat')};
        matlabbatch{1}.spm.stats.ppi.type.ppi.voi = {fullfile(GLM1_dir, ['VOI_' voi '_1.mat'])};
        matlabbatch{1}.spm.stats.ppi.type.ppi.u = [1 1 1; 2 1 -1];
        matlabbatch{1}.spm.stats.ppi.name = PPI_name;
        matlabbatch{1}.spm.stats.ppi.disp = 0;
        % run batch
        disp('Begin to run...');
        spm_jobman('run',matlabbatch);
        disp('PPI specification job completed.');
    end
end


if PPI_GLM_specification == 1
    for ss = loopin
        clear matlabbatch
        
        % specify the current folder
        sub      = num2str(sublist(ss));
        FuncDir  = fullfile(dir_base, sub, func);
        GLM1_dir = fullfile(FuncDir, model, glm1_task);
        
        PPI_dir   = fullfile(GLM1_dir,PPI_folder, PPI_name);
        mkdir(PPI_dir) %output folder
        
        %--------------------------------------------------------------------------------------------------
        clear f files tabfile TAB
        % select functional scans once for all
        f     = spm_select('List', FuncDir, '^s6vol.*\.nii$');
        files = cellstr([repmat([FuncDir '/'],size(f,1),1) f repmat(',1',size(f,1),1)]);
        % select covariate files
        % select experimental condition
        tabfile = fullfile(dir_base, sub,'MNINonLinear/Results/tfMRI_RELATIONAL_LR/RELATIONAL_cleaned_TAB2.txt');
        
        addpath /data/dl577/scripts/spm_scripts
        
        TAB            = importfile_relat(tabfile, 2,37);
        dur_fix        = nanmean(TAB.fix_offset-TAB.fix_onset)/2;
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
        %----------------------------------------------------------------------------------------------------------
        disp(['PPI GLM specification for the No.' num2str(ss) ' subject: ' sub '.'])
        % load the PPI file
        load(fullfile(GLM1_dir,['PPI_' PPI_name '.mat']))
        
        matlabbatch{1}.spm.stats.fmri_spec.dir = {PPI_dir};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 0.72;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        
        
        %%
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = files;
        %%
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
        %
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).name = 'interaction';
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).val = PPI.ppi;
        %
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).name = 'activation';
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val = PPI.P;
        %
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).name = voi;
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).val = PPI.Y;
        %%
            % 1. specify WM.mat files
            load(fullfile(FuncDir, model, 'GLM0', 'VOI_WM_1.mat'))
            WM = Y;
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(4).name = 'WM';
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(4).val = WM;
            clear Y
            % 2. specify CSF.mat files
            load(fullfile(FuncDir, model, 'GLM0', 'VOI_CSF_1.mat'))
            CSF = Y;
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(5).name = 'CSF';
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(5).val = CSF;
            % 3. Block/session effects
            % 3.1
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(6).name = 'BLOCK1';
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(6).val = block1_TRs;
            % 3.2
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(7).name = 'BLOCK2';
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(7).val = block2_TRs;
            % 3.3
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(8).name = 'BLOCK3';
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(8).val = block3_TRs;
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
        disp('PPI GLM specification job completed.');
    end
end

%% Model estimation
if estimation == 1
    for ss = loopin
        clear matlabbatch
        % specify the current folder
        sub      = num2str(sublist(ss));
        FuncDir  = fullfile(dir_base, sub, func);
        GLM1_dir = fullfile(FuncDir, model, glm1_task);
        PPI_dir   = fullfile(GLM1_dir,PPI_folder, PPI_name);
        
        disp(['Estimate PPI glm for the No.' num2str(ss) ' subject: ', sub]);
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {[PPI_dir,'/SPM.mat']}; %Select the SPM.mat of glm1_task
        matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        
        % run batch
        disp('Begin to run...');
        spm_jobman('run',matlabbatch);
        disp('PPI GLM estimation job completed.');
    end
end

%% T Contrast for First-level
if contrast == 1
    for ss = loopin
        % specify the current folder
        sub      = num2str(sublist(ss));
        FuncDir  = fullfile(dir_base, sub, func);
        GLM1_dir = fullfile(FuncDir, model, glm1_task);
        PPI_dir   = fullfile(GLM1_dir,PPI_folder, PPI_name);
        clear matlabbatch
        matlabbatch{1}.spm.stats.con.spmmat = {[PPI_dir,'/SPM.mat']};
        
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = Geffect{1,1};
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        %
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = Geffect{2,1};
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = -1;
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        %
        matlabbatch{1}.spm.stats.con.delete = 1;
        
        % run batch
        disp('Begin to run...');
        disp([Experi ': PPI contrast job completed.']);
        spm_jobman('run',matlabbatch);
    end
end
%toc
%% ------------------------group analysis---------------------- %%

%% GLM2 Specification
if glm2 == 1
    for i = 1:size(Geffect,1)
        gg = Geffect{i,1};
        
        disp(['Group-level specification for the ' gg ' effect.']);
        
        cfiles = {};
        for ss = loopin
            sub      = num2str(sublist(ss));
            FuncDir  = fullfile(dir_base, sub, func);
            GLM1_dir = fullfile(FuncDir, model, glm1_task);
            PPI_dir   = fullfile(GLM1_dir,PPI_folder, PPI_name);
            cfiles{end+1,1} = fullfile( PPI_dir, ['con_000' num2str(Geffect{i,2}) '.nii,1']);
        end
        
        GLM2_dir = fullfile(group_base, model, glm2_task,PPI_folder, PPI_name, gg);
        mkdir(GLM2_dir);
        
        clear matlabbatch
        matlabbatch{1}.spm.stats.factorial_design.dir = {GLM2_dir};
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = cfiles;
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
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
        GLM2_dir = fullfile(group_base, model, glm1_task,PPI_folder, PPI_name, gg);
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
        GLM2_dir = fullfile(group_base, model, glm2_task,PPI_folder, PPI_name, gg);
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
        
        GLM2_dir = fullfile(group_base, model, glm2_task,PPI_folder, PPI_name, gg);
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
