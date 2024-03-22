%% network contructions for 160 ROIs and two hippocampus activations memory similarities
clear all;
clc;
load('leftRestSimilarityLof.mat');
load('rightRestSimilarityLof.mat');
load('ROI160.mat');


Network160={};
allR=[];
for subid = 1:1:24
    ctR=zeros(1,160);
    for sid=1:1:2
        leftHipp=mean(leftRestSimilarity{subid,sid}')';
        rightHipp=mean(rightRestSimilarity{subid,sid}')';
        allTs=[leftHipp,rightHipp,ROI160{subid,sid}];
        R=corrcoef(allTs);
        R(R==1)=0;
        Network160{subid,sid}=R;
        ctR=ctR+R(1,3:1:end)+R(2,3:1:end);
        % allR=[allR;R(2,3:1:end)];
    end
    allR=[allR;ctR./4];
end
% save('Network160.mat','Network160');


[AA,BB,CC,DD]=ttest(allR);
str=spm_vol('rDosenbach160_3mm.nii');
template=spm_read_vols(str);
templateR=zeros(53,63,52);
for i=1:1:160
    xx=find(template==i);
    templateR(xx)=DD.tstat(i);
end
str1=spm_vol('leftTemplate.nii');
str1.fname='ROI160BothHimisphereSimilarityAll.nii';
spm_write_vol(str1,templateR);