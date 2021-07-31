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
func      = 'MNINonLinear/Results/tfMRI_WM_LR/fslsplit'; % functional scan base
glm1_task = 'GLM1_WM_v3';
voi       = {'vPCC', 'dPCC'};
DCM_folder = 'DCM_plus';
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
%}

loopin = 1:3%length(sublist);
dcmmodels = {'fullconn';...
    'VtoD';...
    'DtoV';...
    'VtoD1_DtoV2';...
    'DtoV1_VtoD2';...
    'VtoD1_none2';...
    'DtoV1_none2'};

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
        
        % model 2: vPCC to dPCC in both conditions.
        DCM.b(:,:,1) = [0 0;...  % relat
                        1 0];
        DCM.b(:,:,2) = [0 0; ... % contr
                        1 0];
        
        save(fullfile(DCM_dir,['DCM_' dcmmodels{2} '_' sub '.mat']),'DCM');
        
        % model 3: dPCC to vPCC in both conditions
        clear DCM.b
        DCM.b(:,:,1) = [0 1;... % relat
                        0 0];
        DCM.b(:,:,2) = [0 1;... % contr
                        0 0];
        
        save(fullfile(DCM_dir,['DCM_' dcmmodels{3} '_' sub '.mat']),'DCM');
        
        
        % model 4: vPCC to dPCC in relat, dPCC to vPCC in control
        DCM.b(:,:,1) = [0 0;...
                        1 0];
        DCM.b(:,:,2) = [0 1;...
                        0 0];
        
        save(fullfile(DCM_dir,['DCM_' dcmmodels{4} '_' sub '.mat']),'DCM');
        
        % model 5: dPCC to vPCC in relat, vPCC to dPCC in control
        clear DCM.b
        DCM.b(:,:,1) = [0 1; ...
                        0 0];
        DCM.b(:,:,2) = [0 0; ...
                        1 0];
        
        save(fullfile(DCM_dir,['DCM_' dcmmodels{5} '_' sub '.mat']),'DCM');
        
        % model 6: vPCC to dPCC in relat, none in control
        DCM.b(:,:,1) = [0 0;...
                        1 0];
        DCM.b(:,:,2) = [0 0;...
                        0 0];
        
        save(fullfile(DCM_dir,['DCM_' dcmmodels{6} '_' sub '.mat']),'DCM');
        
        % model 7: dPCC to vPCC in relat, none in control
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
            save(['/data/dl577/matlabbatch/' dcmmodels{m} '_' sub '_fullmodels.mat'], 'sub')
            disp(['DCM estimation (model ' '''' dcmmodels{m} ''')' ' completed for No.' num2str(ss) ' subject: ' sub '.'])
        end       
    end
end

toc

