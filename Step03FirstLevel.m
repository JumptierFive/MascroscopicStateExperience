clear
clc
recycle('on');
rootDir='D:\Work\dataN\ReplayFMRI\';
spm('Defaults','fMRI');
Cnamelist={'a1','a2','a3','a4','a5','b1','b2','b3','b4','b5','c1','c2','c3','c4','c5'};
outputDir=['D:\Work\dataN\ReplayFMRI\FirstLevel'];
for subid=1:1:24
    for SessionId=1:1:2
        for tid=1:1:8
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
                jobs{1}.stats{1}.fmri_spec.dir = {outpath};
                jobs{1}.stats{1}.fmri_spec.timing.units = 'secs';
                jobs{1}.stats{1}.fmri_spec.timing.RT = 2;
                jobs{1}.stats{1}.fmri_spec.timing.fmri_t = 16;
                jobs{1}.stats{1}.fmri_spec.timing.fmri_t0 = 1;
                
                
                
                
                for tid=1
                    jobs{1}.stats{1}.fmri_spec.sess(tid).scans = filename_list(fpath,'PR*.nii');
                    % jobs{1}.stats{1}.fmri_spec.sess(tid).scans = [fpath,'\swratask.nii'];
                    ncon=15;
                    for c=1:ncon
                        jobs{1}.stats{1}.fmri_spec.sess(tid).cond(c).name=Cnamelist{c};
                        jobs{1}.stats{1}.fmri_spec.sess(tid).cond(c).onset=Onset{c}';
                        jobs{1}.stats{1}.fmri_spec.sess(tid).cond(c).duration=3;
                    end
                end
                
                spm_jobman('run',jobs)
                clear jobs
                aaaaaaa
                %******************** Estimate ********************%
                jobs{1}.stats{1}.fmri_est.spmmat = {fullfile(outpath,'SPM.mat')};
                spm_jobman('run',jobs)
                clear jobs
                %******************** Contrast ********************%
                jobs{1}.stats{1}.con.spmmat = {fullfile(outpath,'SPM.mat')};
                jobs{1}.stats{1}.con.consess{1}.fcon.name='effects of interest';
                jobs{1}.stats{1}.con.consess{1}.fcon.convec=eye(16);
                for i=1:15   % 4 is no of conditions
                    jobs{1}.stats{1}.con.consess{i+1}.tcon.name=Cnamelist{i};
                    con=zeros(1,16);con(1,i)=1;
                    jobs{1}.stats{1}.con.consess{i+1}.tcon.convec=con;
                end
                jobs{1}.stats{1}.con.consess{17}.tcon.name='Apha';
                jobs{1}.stats{1}.con.consess{17}.tcon.convec=[1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0];
                jobs{1}.stats{1}.con.consess{18}.tcon.name='Beta';
                jobs{1}.stats{1}.con.consess{18}.tcon.convec=[0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0];
                jobs{1}.stats{1}.con.consess{19}.tcon.name='Gamma';
                jobs{1}.stats{1}.con.consess{19}.tcon.convec=[0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0];
                
                spm_jobman('run',jobs)
                clear jobs
            end
        end
    end
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