clear;
clc;
% load('allSL.mat');

rootDir='D:\Work\dataN\ReplayFMRI\FirstLevel';
allFiles=filename_list(rootDir,'sub*');
output=['D:\Work\dataN\ReplayFMRI\SimilarityForAll'];mkdir(output);
tempN=spm_read_vols(spm_vol('template.nii'));
tempN=reshape(tempN,[1,53*63*52]);

for i=1:1:length(allFiles)
    for si=1:1:2
        
        outputDir=[output,'\subid',dec2base(i,10,2),'\session',dec2base(si,10,2)];mkdir(outputDir);
        for mid=1:1:15
            template=[];
            
            for taskid=1:1:8
                templateH=spm_vol([allFiles{i},'\Session',dec2base(si,10,2),'\task_',dec2base(taskid,10,2),'\taskid_',dec2base(mid,10,3),'.nii']);
                templateN=spm_read_vols(templateH);
                template=[template;reshape(templateN,[1,53*63*52])];
            end
            [~,~,~,DD]=ttest(template);
            templateH.fname=[outputDir,'\memid',dec2base(mid,10,3),'.nii'];
            tValue=DD.tstat;tValueR=tValue;tValueR(isnan(tValueR))=[];tValueR(tValueR<-100)=-100;tValueR(tValueR>100)=100;
            meanT=mean(tValueR);sdT=std(tValueR);maxT=meanT+sdT*5;minT=meanT-sdT*5;
            maxT
            minT
            tValue(tValue>maxT)=maxT;tValue(tValue<minT)=minT;
            spm_write_vol(templateH,reshape(tValue,[53,63,52]));
        end
        
        
    end
end