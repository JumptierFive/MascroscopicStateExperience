clear all;
clc;
nii_info = spm_vol('AAL116_1mm.nii');
original_volume = spm_read_vols(nii_info);

% ����ԭʼ��Ŀ������سߴ�
original_size = size(original_volume);
target_size = [53, 63, 52];

% �������ű���
scale_factors = target_size ./ original_size;

% ʹ��imresize���в�ֵ
resized_volume = imresize3(original_volume, target_size, "nearest");

fileN=spm_vol('template.nii');
fileN.fname='resized.nii';
spm_write_vol(fileN,resized_volume);