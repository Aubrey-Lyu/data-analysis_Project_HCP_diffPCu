%% %%%%%%%%% HCP wm/n-back task ROI timeseries extraction respectively for dPCu and vPCu %%%%%%%%%%%%
%-----------------------------------------------------------------------
clear
%% What subjects?
%Experi = 'HCP_wm_activation';
sublist = [100307
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

%% What analysis
contrast = 1;
VOI_extr = 1;
%% Define directories
masks    = {'vPCC','vPCC.nii';'dPCC','dPCC.nii'};
MaskDir   = '....';
dir_base  = '...';
func      = '....'; % functional scan base
model     = 'model';
glm1_task = 'GLM1_WM_v3';

loopin = 1:length(sublist);

%Initiating SPM
%spm fmri
spm_jobman('initcfg');
%% Contrast to create F contrast
if contrast == 1
    for ss = loopin
        % specify the current folder
        sub      = num2str(sublist(ss));
        FuncDir  = fullfile(dir_base, sub, func);
        GLM1_dir   = fullfile(FuncDir, model, glm1_task);
        
        F = [];
        % how to define F contrast matrix
        load(fullfile(GLM1_dir, 'SPM.mat'))
        maineffect = length(SPM.Sess.U);
        alleffect = length(SPM.Sess.col);
        %num_con = length(SPM.xCon);
        
        F = [eye(maineffect), zeros(maineffect, alleffect-maineffect)]; % adjust for main effects - F contrast
        clear matlabbatch
        
        matlabbatch{1}.spm.stats.con.spmmat = {[GLM1_dir,'/SPM.mat']};
        matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'main effects';
        matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights = F;
        matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.delete = 0;
        % run batch
        disp('Begin to run...');
        spm_jobman('run',matlabbatch);
        disp('glm1 contrast job completed.');
        
    end
end


%% VOIs Extraction

if VOI_extr ==1
    for ss = loopin
        % specify the current folder
        sub      = num2str(sublist(ss));
        FuncDir  = fullfile(dir_base, sub, func);
        GLM1_dir   = fullfile(FuncDir, model, glm1_task);
        for m = 1:size(masks,1)
            
            clear matlabbatch SPM
            
            maskname = masks{m,1};
            mask     = masks{m,2};
            disp(['Extracting timeseries of ' maskname ' for No.' num2str(ss) ' subject: ' sub '.'])
            load(fullfile(GLM1_dir, 'SPM.mat'))
            
            % find the index of the F contrast
            for i = 1:length(SPM.xCon)
                if strcmp(SPM.xCon(i).STAT, 'F')
                    F_ind = i;
                    break
                end
            end
            
            matlabbatch{1}.spm.util.voi.spmmat = {[GLM1_dir, '/SPM.mat']};
            matlabbatch{1}.spm.util.voi.adjust = F_ind;
            matlabbatch{1}.spm.util.voi.session = 1;
            matlabbatch{1}.spm.util.voi.name = maskname;
            matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {fullfile(MaskDir, mask)};
            matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.5;
            matlabbatch{1}.spm.util.voi.expression = 'i1';
            
            % run batch
            disp('Begin to run...');
            spm_jobman('run',matlabbatch);
            
        end
    end
    if ss == length(sublist)
        disp('Job completed for the group.');
    end
end
