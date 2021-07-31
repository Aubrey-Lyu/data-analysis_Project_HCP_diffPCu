%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DYNAMIC CAUSAL MODELLING for task_HCP relational task
% for the second kind of model: cross influence for relat vs contr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
rmpath(genpath('/applications/spm/spm12_6906'))
rmpath(genpath('/applications/spm/spm12_7219'))
addpath('/home/dl577/spm12')
addpath('/data/dl577/scripts/spm_scripts/')

%% What analysis
dcm_spec = 0;
dcm_esti = 1;
%% Define directories
dir_base  = '/lustre/scratch/wbic-beta/dl577/HCP100/';
group_base= '/data/dl577/task_HCP/HCP100/group_results';
model     = 'model';
func      = 'MNINonLinear/Results/tfMRI_RELATIONAL_LR/fslsplit_scans'; % functional scan base
glm1_task = 'GLM1_RELATIONAL_v3';
voi       = {'aPCC', 'pPCC'};
DCM_folder = 'DCM_plus';
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
%}

loopin = 1:15%length(sublist);
dcmmodels = {'fullconn';...
    'AtoP';...
    'PtoA';...
    'AtoP1_PtoA2';...
    'PtoA1_AtoP2';...
    'AtoP1_none2';...
    'PtoA1_none2'};

%% Loop through subjects

if dcm_spec == 1
    for ss = loopin
        clear matlabbatch
        % specify the current folder
        sub      = num2str(sublist(ss));
        FuncDir  = fullfile(dir_base, sub, func);
        GLM1_dir = fullfile(FuncDir, model, glm1_task);
        
        
        disp(['Specify DCM Model for No.' num2str(ss) ' subject: ', sub '.']);
        
        %load(fullfile(GLM1_dir,'SPM.mat'));
        DCM_dir = fullfile(GLM1_dir, DCM_folder);
        mkdir(DCM_dir)
        %
        %load('/data/dl577/task_HCP/Temp_Variables/metaB.mat')
        load(fullfile(GLM1_dir,'SPM.mat'))
        % Load regions of interest
        %--------------------------------------------------------------------------
        for i = 1:length(voi)
            clear xY
            load(fullfile(GLM1_dir,['VOI_' voi{i} '_1.mat']),'xY');
            DCM.xY(i) = xY;
        end
        clear i
        
        DCM.n = length(DCM.xY);      % n = number of regions
        DCM.v = length(DCM.xY(1).u); % v = number of time points
        
        % Time series
        %--------------------------------------------------------------------------
        DCM.Y.dt  = SPM.xY.RT;
        DCM.Y.X0  = DCM.xY(1).X0;
        for i = 1:DCM.n
            DCM.Y.y(:,i)  = DCM.xY(i).u;
            DCM.Y.name{i} = DCM.xY(i).name;
        end
        
        DCM.Y.Q    = spm_Ce(ones(1,DCM.n)*DCM.v);
        
        % Experimental inputs
        %--------------------------------------------------------------------------
        DCM.U.dt   =  SPM.Sess(1).U(1).dt;
        DCM.U.name = [SPM.Sess.U.name];
        DCM.U.u    = [SPM.Sess(1).U(1).u(33:end,1) ...
            SPM.Sess(1).U(2).u(33:end,1)];
        
        % DCM parameters and options
        %--------------------------------------------------------------------------
        DCM.delays = repmat(SPM.xY.RT/2,DCM.n,1);
        DCM.TE     = 0.0331;
        
        DCM.options.nonlinear  = 1;
        DCM.options.two_state  = 1;
        DCM.options.stochastic = 1;
        DCM.options.nograph    = 1;
        
        %
        % Connectivity matrices for model with backward modulation
        %--------------------------------------------------------------------------
        DCM.a = ones(2,2); % Intrinsic connectivity: fully connected
        DCM.c = zeros(2,2);  % Direct input to specified regions: no
        DCM.d = zeros(2,2,0); % Modulatory effect on specified connectivity: no
        
        %%  B matrix specification
        % model 1: fully connected
        DCM.b(:,:,1) = [0 1;... % relat
                        1 0];
        DCM.b(:,:,2) = [0 1;... % contr
                        1 0];
        
        save(fullfile(DCM_dir,['DCM_' dcmmodels{1} '_' sub '.mat']),'DCM');
        
        % model 2: aPCC to pPCC in both conditions.
        DCM.b(:,:,1) = [0 0;...  % relat
                        1 0];
        DCM.b(:,:,2) = [0 0; ... % contr
                        1 0];
        
        save(fullfile(DCM_dir,['DCM_' dcmmodels{2} '_' sub '.mat']),'DCM');
        
        % model 3: pPCC to aPCC in both conditions
        clear DCM.b
        DCM.b(:,:,1) = [0 1;... % relat
                        0 0];
        DCM.b(:,:,2) = [0 1;... % contr
                        0 0];
        
        save(fullfile(DCM_dir,['DCM_' dcmmodels{3} '_' sub '.mat']),'DCM');
        
        
        % model 4: aPCC to pPCC in relat, pPCC to aPCC in control
        DCM.b(:,:,1) = [0 0;...
                        1 0];
        DCM.b(:,:,2) = [0 1;...
                        0 0];
        
        save(fullfile(DCM_dir,['DCM_' dcmmodels{4} '_' sub '.mat']),'DCM');
        
        % model 5: pPCC to aPCC in relat, aPCC to pPCC in control
        clear DCM.b
        DCM.b(:,:,1) = [0 1; ...
                        0 0];
        DCM.b(:,:,2) = [0 0; ...
                        1 0];
        
        save(fullfile(DCM_dir,['DCM_' dcmmodels{5} '_' sub '.mat']),'DCM');
        
        % model 6: aPCC to pPCC in relat, none in control
        DCM.b(:,:,1) = [0 0;...
                        1 0];
        DCM.b(:,:,2) = [0 0;...
                        0 0];
        
        save(fullfile(DCM_dir,['DCM_' dcmmodels{6} '_' sub '.mat']),'DCM');
        
        % model 7: pPCC to aPCC in relat, none in control
        clear DCM.b
        DCM.b(:,:,1) = [0 1; ...
                        0 0];
        DCM.b(:,:,2) = [0 0; ...
                        0 0];
        
        save(fullfile(DCM_dir,['DCM_' dcmmodels{7} '_' sub '.mat']),'DCM');
        
        disp(['DCMs specification completed for No.' num2str(ss) ' subject: ' sub '.'])
    end
end

tic
if dcm_esti == 1
    for ss = loopin
        % specify the current folder
        sub      = num2str(sublist(ss));
        FuncDir  = fullfile(dir_base, sub, func);
        GLM1_dir = fullfile(FuncDir, model, glm1_task);
        DCM_dir = fullfile(GLM1_dir, DCM_folder);
        
        for m = 1:length(dcmmodels)
            spm_dcm_estimate(fullfile(DCM_dir,['DCM_' dcmmodels{m} '_' sub '.mat']));
            save(['/data/dl577/matlabbatch/' dcmmodels{m} '_' sub '.mat'], 'sub')
            disp(['DCM estimation (model ' '''' dcmmodels{m} ''')' ' completed for No.' num2str(ss) ' subject: ' sub '.'])
        end       
    end
end

toc

