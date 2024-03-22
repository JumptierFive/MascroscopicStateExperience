clear all;
clc;
str=spm_vol('template.nii');
template=spm_read_vols(str);


for x=1:1:27
    for y=1:1:63
        for z=1:1:52
            template(x,y,z)=0;
        end
    end
end

str.fname='LeftTemplate.nii';
spm_write_vol(str,template);