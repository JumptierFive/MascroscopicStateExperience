%% preprocessing for the functional activations
clear all;
clc;
rootDir='D:\Work\dataN\ReplayFMRI';
for subid=1:1:24
    for sid=1:1:2
        for tid=1:1:8
            fpath=[rootDir,'\subid',dec2base(subid,10,2),'\Session',dec2base(sid,10,2),'\task_',dec2base(tid,10,2)];
            cd(fpath);
            %% slice timing
            fileNamelist=filename_list(fpath,'*.nii');
            %% fileNamelist=fileNamelist(11:end);
            jobs{1}.temporal{1}.st.scans = {fileNamelist};
            jobs{1}.temporal{1}.st.nslices = 36;
            jobs{1}.temporal{1}.st.tr = 2;
            jobs{1}.temporal{1}.st.ta =2-2/36;
            jobs{1}.temporal{1}.st.so = [2:2:36 1:2:35];
            jobs{1}.temporal{1}.st.refslice = 1;
            jobs{1}.temporal{1}.st.prefix = 'a';
            spm_jobman('run',jobs)
            clear jobs
            
            %% head motion correction
            jobs{1}.spatial{1}.realign{1}.estwrite.data = {filename_list(fpath,'a*.nii') };
            jobs{1}.spatial{1}.realign{1}.estwrite.roptions.prefix = 'r';
            spm_jobman('run',jobs)
            clear jobs
            
            
            
            %% segment
            jobs{1}.spatial{1}.preproc.channel.vols = filename_list(fpath,'meana*.nii');
            jobs{1}.spatial{1}.preproc.warp.affreg = 'mni';
            jobs{1}.spatial{1}.preproc.warp.write = [0 1];
            spm_jobman('run',jobs)
            clear jobs
            
            %% normalize
            % jobs{1}.spatial{1}.normalise{1}.write.subj.def = filename_list(spath,'y*.nii');
            jobs{1}.spatial{1}.normalise{1}.write.subj.def = filename_list(fpath,'y*.nii');
            jobs{1}.spatial{1}.normalise{1}.write.subj.resample = filename_list(fpath,'r*.nii');
            jobs{1}.spatial{1}.normalise{1}.write.woptions.bb = [-90 -126 -72 90 90 108];   % two point in MNI
            jobs{1}.spatial{1}.normalise{1}.write.woptions.vox = [3 3 3];   % voxel volume
            jobs{1}.spatial{1}.normalise{1}.write.woptions.interp = 1;
            spm_jobman('run',jobs)
            clear jobs
            
            
            % spatial moothing
            jobs{1}.spatial{1}.smooth.data = filename_list(fpath,'w*.nii');
            jobs{1}.spatial{1}.smooth.fwhm = [6 6 6];
            spm_jobman('run',jobs)
            clear jobs
        end
    end
end