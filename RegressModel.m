%% 
clear all;
clc;
allInfo=readtable('Memo.xls');
y=allInfo.Memory1Verbal;
X=[ones(24,1),allInfo.age,allInfo.sexN];
[b,bint,r1,rint,stats] = regress(y,X);

y=allInfo.Memory1Visual;
X=[ones(24,1),allInfo.age,allInfo.sexN];
[b,bint,r2,rint,stats] = regress(y,X);

rr=r1+r2;