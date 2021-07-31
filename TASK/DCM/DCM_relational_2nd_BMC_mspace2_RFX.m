%-----------------------------------------------------------------------
% Job saved on 28-Oct-2018 13:37:14 by dl577
% DCM - second level - Baysian Model Comparison
% Data: HCP - relational task - 100 subjs
%-----------------------------------------------------------------------
clear

sublist=[188347
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

%% Change directories here
base_dir  = '/home/dl577/scratch/HCP100/';
Funcdir   = '/MNINonLinear/Results/tfMRI_RELATIONAL_LR/fslsplit_scans/model/GLM1_RELATIONAL_v3/DCM_plus/';
outputdir = '/data/dl577/task_HCP/HCP100/DCM/model_space2_RFX/relational';
mkdir(outputdir)
dcmmodels = {'fullconn';...
    'AtoP';...
    'PtoA';...
    'AtoP1_PtoA2';...
    'PtoA1_AtoP2';...
    'AtoP1_none2';...
    'PtoA1_none2'}; % 'AtoP': vPCu --> dPCu for both conditions; 'AtoP1_PtoA2': vPCu->dPCu in difficult condition, dPCu->vPCu in easy condition

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


