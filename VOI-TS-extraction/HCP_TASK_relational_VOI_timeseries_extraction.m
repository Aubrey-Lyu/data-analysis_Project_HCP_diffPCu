%% %%%%%%%%% HCP relational task ROI timeseries extraction respectively for dPCu and vPCu %%%%%%%%%%%%
%-----------------------------------------------------------------------
clear
%% What subjects?
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
%% What analysis
contrast = 1;
VOI_extr = 1;
%% Define directories
masks    = {'vPCC','vPCC.nii';'dPCC','dPCC.nii'}; % dPCu, vPCu as reported names
MaskDir   = '...';
dir_base  = '...';
func      = '...'; % functional scan base
model     = 'model';
glm1_task = 'GLM1_RELATIONAL_v3';

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
