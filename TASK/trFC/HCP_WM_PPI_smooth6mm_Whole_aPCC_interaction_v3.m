clear
rmpath(genpath('/applications/spm/spm12_6906'))
rmpath(genpath('/applications/spm/spm12_7219'))
addpath('/home/dl577/spm12')
addpath('/data/dl577/scripts/spm_scripts/')
Experi = 'taskWM_PPI_6mm_aPCC_interaction_';
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
sublist  = [100307
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
    899885];

loopin = 2:length(sublist);

%% Define directories
dir_base  = '/lustre/scratch/wbic-beta/dl577/HCP100/';
group_base= '/data/dl577/task_HCP/HCP100/group_results';
model     = 'model';
func      = 'MNINonLinear/Results/tfMRI_WM_LR/fslsplit'; % functional scan base
struc     = 'T1w/Results/tfMRI_WM_LR';  % structural scan base
glm1_task = 'GLM1_WM_v3';
glm2_task = 'GLM2_WM_v3';
voi       = 'vPCC';
PPI_folder = 'PPI_analyses';
PPI_name  = [voi '_interaction'];
Geffect   = {[voi '_2bk-0bk'],1;...
            [voi '_0bk-2bk'], 2};

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
        TAB = importfile(fullfile(dir_base, sub,'MNINonLinear/Results/tfMRI_WM_LR/WM_run2_TAB.txt'));

        % select session break point
        fix_onset = (TAB.Fix15secOnsetTime(strcmp(TAB.ProcedureBlock,'Fix15secPROC'))-TAB.SyncSlideOnsetTime(1))/1000;
        fix_ind = find(strcmp(TAB.ProcedureBlock,'Fix15secPROC'));
        cue_ind = fix_ind + 1;
        events = [TAB.Cue2BackOnsetTime; TAB.CueTargetOnsetTime; TAB.Fix15secOnsetTime; TAB.StimOnsetTime];
        events = ([NaN(4,1); sort(events(~isnan(events)))]-TAB.SyncSlideOnsetTime(1))/1000; % the first 4 NaN correspond to the dummy scans, so to sync with the TAB.ProcedureBlock column
        dur = (events(cue_ind(1:end-1)) - events(fix_ind(1:end-1)))/2;
        session_breaks = events(fix_ind(1:end-1)) + dur;
        tr_breaks     = round(session_breaks/0.72, 0);
        final_break   = round((fix_onset(4)+mean(dur))/0.72, 0);
        trs            = zeros(length(files),1);
        block1_TRs     = trs;
        block1_TRs(5:tr_breaks(1)) = 1;
        block2_TRs     = trs;
        block2_TRs(tr_breaks(1):tr_breaks(2)) = 1;
        block3_TRs     = trs;
        block3_TRs(tr_breaks(2):final_break) = 1;
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
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {fullfile(dir_base,sub,'MNINonLinear/Results/tfMRI_WM_LR/Movement_Regressors.txt')};
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
        for ss = 1:length(sublist)
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
        GLM2_dir = fullfile(group_base, model, glm2_task,PPI_folder, PPI_name, gg);
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
