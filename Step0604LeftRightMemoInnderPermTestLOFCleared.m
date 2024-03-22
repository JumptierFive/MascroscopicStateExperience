%% compare activation similarity between real and dummy sitations
clear all;
clc;
load('RightHippVols.mat');
HippR=allHippVols;
load('RightHippVolsPerm.mat');
allSimi={};
for subid=1:1:24
    subid
    for sid=1:1:2
        allX=[];
        for mid=1:1:15
            mid
            xn=0;
            nn=0;
            for tidx=1:1:8
                for tidy=1:1:8
                    if tidx==tidy
                        continue;
                    end
                    R=corrcoef(HippR{subid,sid}{mid}(tidx,:)',HippR{subid,sid}{mid}(tidy,:)');
                    R=R(2);
                    Distri=[];
                    for mii=1:1:15
                        for miii=1:1:15
                            T=corrcoef(HippR{subid,sid}{mii}(tidx,:)',HippR{subid,sid}{miii}(tidy,:)');
                            Distri=[Distri;T(2)];
                        end
                    end
                    Distri(isnan(Distri))=[];
                    xx=find(Distri<R);
                    xx=numel(xx)/numel(Distri);
                    xn=xn+xx;
                    nn=nn+1;
                end
            end
            
            allX=[allX;xn/nn];
        end
        allSimi{subid,sid}=allX;
    end
end

rightSubSimi=[];
for subid=1:1:24
    a1=[mean(allSimi{subid,1}),mean(allSimi{subid,2})];
    rightSubSimi=[rightSubSimi;a1];
end



load('LeftHippVols.mat');
HippR=allHippVols;
load('LeftHippVolsPerm.mat');
allSimi={};
for subid=1:1:24
    subid
    for sid=1:1:2
        allX=[];
        for mid=1:1:15
            mid
            xn=0;
            nn=0;
            for tidx=1:1:8
                for tidy=1:1:8
                    if tidx==tidy
                        continue;
                    end
                    R=corrcoef(HippR{subid,sid}{mid}(tidx,:)',HippR{subid,sid}{mid}(tidy,:)');
                    R=R(2);
                    Distri=[];
                    for mii=1:1:15
                        for miii=1:1:15
                            T=corrcoef(HippR{subid,sid}{mii}(tidx,:)',HippR{subid,sid}{miii}(tidy,:)');
                            Distri=[Distri;T(2)];
                        end
                    end
                    Distri(isnan(Distri))=[];
                    xx=find(Distri<R);
                    xx=numel(xx)/numel(Distri);
                    xn=xn+xx;
                    nn=nn+1;
                end
            end
            
            allX=[allX;xn/nn];
        end
        allSimi{subid,sid}=allX;
    end
end

leftSubSimi=[];
for subid=1:1:24
    a1=[mean(allSimi{subid,1}),mean(allSimi{subid,2})];
    leftSubSimi=[rightSubSimi;a1];
end