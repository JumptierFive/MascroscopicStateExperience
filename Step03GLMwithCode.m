clear
clc
recycle('on');
rootDir='D:\dataN\ReplayFMRI\';
spm('Defaults','fMRI');
Cnamelist={'a1','a2','a3','a4','a5','b1','b2','b3','b4','b5','c1','c2','c3','c4','c5'};
outputDir=['D:\dataN\ReplayFMRI\FirstLevel'];
hrf = spm_hrf(1);
headFF=spm_vol('D:\dataN\Session01\sub01RUN01\bcWGSdswranrun01.nii');
headFF=headFF(1,:);
for subid=1:1:24
    subid
    for SessionId=1:1:2
        SessionId
        for tid=1:1:8
            tid
            fpath=[rootDir,'\Processed\subid',dec2base(subid,10,2),'\Session',dec2base(SessionId,10,2),'\task_',dec2base(tid,10,2),'\'];
            outpath=[outputDir,'\subid',dec2base(subid,10,2),'\Session',dec2base(SessionId,10,2),'\task_',dec2base(tid,10,2),'\'];
            mkdir(outpath);
            
            if exist(fpath,'dir')
                [currentH,currentN,currentTxt] = tsvread([fpath,'/event.tsv']);
                currentH(1,:)=[];currentH=currentH(:,1);currentTxt(1,:)=[];currentTxt=currentTxt(:,3);
                Onset={};
                for mi=1:1:15
                    Onset{mi}=[];
                end
                for oi=1:1:length(currentH)
                    if strcmp(currentTxt{oi},'alpha1')
                        Onset{1}=[Onset{1};currentH(oi)];
                    end
                    if strcmp(currentTxt{oi},'alpha2')
                        Onset{2}=[Onset{2};currentH(oi)];
                    end
                    if strcmp(currentTxt{oi},'alpha3')
                        Onset{3}=[Onset{3};currentH(oi)];
                    end
                    if strcmp(currentTxt{oi},'alpha4')
                        Onset{4}=[Onset{4};currentH(oi)];
                    end
                    if strcmp(currentTxt{oi},'alpha5')
                        Onset{5}=[Onset{5};currentH(oi)];
                    end
                    
                    
                    if strcmp(currentTxt{oi},'beta1')
                        Onset{6}=[Onset{6};currentH(oi)];
                    end
                    if strcmp(currentTxt{oi},'beta2')
                        Onset{7}=[Onset{7};currentH(oi)];
                    end
                    if strcmp(currentTxt{oi},'beta3')
                        Onset{8}=[Onset{8};currentH(oi)];
                    end
                    if strcmp(currentTxt{oi},'beta4')
                        Onset{9}=[Onset{9};currentH(oi)];
                    end
                    if strcmp(currentTxt{oi},'beta5')
                        Onset{10}=[Onset{10};currentH(oi)];
                    end
                    
                    
                    if strcmp(currentTxt{oi},'gamma1')
                        Onset{11}=[Onset{11};currentH(oi)];
                    end
                    if strcmp(currentTxt{oi},'gamma2')
                        Onset{12}=[Onset{12};currentH(oi)];
                    end
                    if strcmp(currentTxt{oi},'gamma3')
                        Onset{13}=[Onset{13};currentH(oi)];
                    end
                    if strcmp(currentTxt{oi},'gamma4')
                        Onset{14}=[Onset{14};currentH(oi)];
                    end
                    if strcmp(currentTxt{oi},'gamma5')
                        Onset{15}=[Onset{15};currentH(oi)];
                    end
                end
                
                fileName=[fpath,'\swratask.nii'];
                volsM=spm_vol(fileName);
                volsN=spm_read_vols(volsM);
                
                timeSO={};
                templateAll={};
                for mid=1:1:15
                    onsetCt=Onset{mid}./2;
                    timeS=zeros(length(volsM),1);
                    for tii=1:1:length(timeS)
                        ann=tii-onsetCt;
                        ann(ann<0)=[];
                        if min(ann)<=3
                            timeS(tii)=1;
                        end
                    end
                    
                    timeSN=conv(timeS,hrf);
                    timeSN=timeSN(1:1:length(timeS));
                    timeSO{mid}=timeSN;
                    template=zeros(53,63,52);
                    templateAll{mid}=template;
                end
                
                
                for zi=1:1:53
                    for zzi=1:1:63
                        for zzii=1:1:52
                            AD=volsN(zi,zzi,zzii,:);
                            AD=AD(:);
                            if sum(AD)==0||std(AD)==0 || numel(unique(AD))<=5
                                continue;
                            end
                            for mmid=1:1:15
                                % beta=corrcoef(AD,timeSO{mmid});
                                beta=AD\timeSO{mmid};
                                resd=AD-timeSO{mmid}*beta;
                                t_values=beta./std(resd);
                                templateAll{mmid}(zi,zzi,zzii)=t_values;
                            end
                        end
                    end
                end
                for miid=1:1:15
                    headFF.fname=[outpath,'\taskid_',dec2base(miid,10,3),'.nii'];
                    spm_write_vol(headFF,templateAll{miid});
                end
                
            end
        end
    end
    clc
end





function varargout = tsvread( varargin )
%[data, header, raw] = tsvread( file ) reads in text file with tab-seperated variables. default value for data is nan.
%alternative input/output option is suppluying header strings
%[col1, col2, col3, ..., header, raw] = tsvread( file, header1, header2, header3, ... )
%header is the first row (assumed to have header names) and raw is the imported text
%if a vector is supplied, this specifies the number rows to be imported.
%examples:
%[col1, col2,col3,header,raw] = tsvread( 'example.tsv', 'header1', 'header2', 'header3', 1:5 )
%will import data from example.tsv, and cols corresponding to header1,
%header2, header3, cols 1 t0 5.
%if no outputs are requested, then a portion of the rquested table is
%displayed -- good idea to see how the import and header requests are
%working!
%If there is a tsv file in local directory, then just running "tsvread" at
%the command line will read the newest tsv file and display first ten rows to screen.
%sak 9/2/11
if nargin == 0
    disp( 'importing most recent tsv file in local directory' );
    d = dir( '*.tsv' );
    [~,i] = sort( [d.datenum], 'descend' );
    fprintf( 'tsvread( ''%s'', 1:10 )\n', d(i(1)).name );
    eval( sprintf( 'tsvread( ''%s'', 1:10 )', d(i(1)).name ) );
    return;
end;
fid = fopen( varargin{1}, 'r' );
if nargout == 0
    fprintf( 'loading %s, and displaying sideways\n', varargin{1} );
end
varargin(1) = [];
stuff = textscan( fid, '%s', 'delimiter', '\n');
stuff = stuff{1};
fclose(fid);
numrows = 1:size(stuff,1);
ind = cellfun( 'isclass', varargin, 'double' );
if any( ind )
    numrows = varargin{ind};
    varargin(ind) = [];
    if numel(numrows) == 1
        numrows = 1:min(numrows, size( stuff, 1) );
    end
end
numrows = intersect( numrows, 1:size( stuff, 1 ) );
header = regexprep( regexp( stuff{1}, '[^\t]*\t', 'match' ), '\t', '' );
raw = repmat( {}, numel(numrows), numel(header) );
for i=numrows
    stuff{i}(end+1) = 9;
    tmp = regexprep( regexp( stuff{i}, '[^\t]*\t', 'match' ), '\t', '');
    raw(i,1:numel(tmp)) = tmp;
end
header = raw(1,:);
data = nan(size(raw));
for i=numrows
    for j=1:size( raw, 2 )
        if ~isempty( raw{i,j} )
            [a, count, errmsg] = sscanf( raw{i,j}, '%f' );
            if ~isempty( a )
                data(i,j) = a;
            end
        end
    end
end
%%
if numel( varargin ) == 0
    j=1:size(data, 2);
else
    j = [];
    for i=1:numel(varargin)
        j = [j, find( strncmp( header, varargin{i}, numel(varargin{i})) )];
    end
end
if nargout == 0
    disp( [ strvcat( header(j)), num2str( data(numrows,j)' ) ] );
    return;
end
if numel(varargin)==0
    varargout = {data, header, raw};
    varargout = varargout(1:nargout);
    return;
end
varargout = {};
for i=j
    varargout{end+1} = data( :, i );
end
varargout{end+1} = header(:,j);
varargout{end+1} = raw(:,j);
end