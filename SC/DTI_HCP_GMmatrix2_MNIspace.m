clear
% convert the streamline count (terminating in Grey matter) stored in the matrix format to the brain space/nifti format

base_dir = '/lustre/scratch/wbic-beta/dl577/HCP100/DTI_analysis';

for g = 3:9
    group_folder = ['group' num2str(g)];
    for voi_folder = {'dPre_DiffSpace', 'vPre_DiffSpace'}
        voi_folder = char(voi_folder);
        sub_folders = dir(fullfile(base_dir, group_folder));
        for s = 3:length(sub_folders) % count from the third to the end
            sub_folder = sub_folders(s).name;
          
            data_dir =  fullfile(base_dir, group_folder, sub_folder, 'T1w/Diffusion.bedpostX', voi_folder);  
            disp(['Subj ' sub_folder '-' voi_folder ': Preparing to write image...'])
            % Load Matrix2   
            clear x M GM_mean_track coor data ind
            x=load(fullfile(data_dir, 'fdt_matrix2.dot'));
            M=full(spconvert(x));
            GM_mean_track = mean(M);
            coor=load(fullfile(data_dir, 'tract_space_coords_for_fdt_matrix2'))+1; % +1 to convert fsl naming convention
            % Create a nifi file with the data filled in
            data = zeros(145,174,145);
            % using vector indices is a little faster than using loop
            ind = sub2ind(size(data), coor(:,1), coor(:,2), coor(:,3));
            data(ind) = GM_mean_track;
            
            %% Write nifti to a standard img
            addpath('~/scripts/spm_scripts/NIfTI_20140122/NIfTI_20140122')
            standard_img = '/lustre/scratch/wbic-beta/dl577/HCP100/DTI_analysis/group0/100307/T1w/T1w_acpc_dc_restore_1.25.nii.gz';
            [img_matrix, hdr] = generic_script(standard_img);
            % Create nifti file data structure
            nii.img = data;
            nii.hdr = hdr;
            save_nii(nii, fullfile(data_dir,'GM_mean_matrix2.nii'))%
            disp(['Subj ' sub_folder '-' voi_folder '-- Image is written.'])
           
        end
    end
end
