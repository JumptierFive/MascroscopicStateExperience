%% replay step 01, data files collection for GLM files
%% Leinian Li 2023 12 15 at Jinan, Shandong Normal Univerisity
clear all;clc;
rootDir='D:\Work\dataN';
allFiles=filename_list([rootDir,'\rawData'],'sub-*');

for i=1:1:length(allFiles)
    ctFile=filename_list([allFiles{i},'\ses-2\func'],'*task-satellite_run-0*bold.nii.gz');
    ctFileB=filename_list([allFiles{i},'\ses-2\func'],'*task-satellite_run-0*_events.tsv');
    for ki=1:1:length(ctFile)
        outputDir=[rootDir,'\ReplayFMRI\','\subid',dec2base(i,10,2),'\Session02\task_',dec2base(ki,10,2),'\'];mkdir(outputDir);
        targetFile=[outputDir,'\task.nii.gz'];
        copyfile(ctFile{ki},targetFile);
        gunzip(targetFile);
        targetFile2=[outputDir,'\event.tsv'];
        copyfile(ctFileB{ki},targetFile2);
    end
end