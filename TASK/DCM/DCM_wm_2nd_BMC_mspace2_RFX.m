%-----------------------------------------------------------------------
% Job saved on 28-Oct-2018 13:37:14 by dl577
% DCM - second level - Baysian Model Comparison
% Data: HCP - relational task - 100 subjs
%-----------------------------------------------------------------------
clear

%-----------------------------------------------------------------------
% Job saved on 28-Oct-2018 13:37:14 by dl577
% DCM - second level - Baysian Model Comparison
% Data: HCP - WM task - 100 subjs
%-----------------------------------------------------------------------
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
loopin = 1:length(sublist);

%% directories
base_dir = '/home/dl577/scratch/HCP100/';
Funcdir = '/MNINonLinear/Results/tfMRI_WM_LR/fslsplit/model/GLM1_WM_v3/DCM_plus/';
DCM_folder = {'DCM', 'DCM_plus'};
outputdir = '/data/dl577/task_HCP/HCP100/DCM/model_space2_RFX/wm/';
mkdir(outputdir)

dcmmodels = {'fullconn';...
    'VtoD';...
    'DtoV';...
    'VtoD1_DtoV2';...
    'DtoV1_VtoD2';...
    'VtoD1_none2';...
    'DtoV1_none2'};

%% ----------- matlabbatch ---------------
clear matlabbatch

matlabbatch{1}.spm.dcm.bms.inference.dir = {outputdir};

for ss = loopin
    sub = num2str(sublist(ss));
    models = {}; % a cell to feed in matlabbatch, directing to all the DCM files of ONE subject.
    for m = 1:length(dcmmodels)
        dcmmodel = dcmmodels{m};
        models{end+1,1} = [base_dir sub Funcdir 'DCM_' dcmmodel '_' sub '.mat'];
    end
matlabbatch{1}.spm.dcm.bms.inference.sess_dcm{ss}.dcmmat = models;
end

matlabbatch{1}.spm.dcm.bms.inference.model_sp = {''};
matlabbatch{1}.spm.dcm.bms.inference.load_f = {''};
matlabbatch{1}.spm.dcm.bms.inference.method = 'RFX';
matlabbatch{1}.spm.dcm.bms.inference.family_level.family_file = {''};
matlabbatch{1}.spm.dcm.bms.inference.bma.bma_yes.bma_famwin = 'famwin';
matlabbatch{1}.spm.dcm.bms.inference.verify_id = 1;

% run batch
disp('Begin to run...');
spm_jobman('run',matlabbatch);

