%https://sites.google.com/site/satayusers/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Supplementary script 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all

[file, path]=uigetfile('*.bam'); %<-point to the BAM file
cd(path)

%%%load some variables
infobam=baminfo(file,'ScanDictionary',true);
load('yeastGFF.mat') %<- both yeastGFF.mat should be in the same folder as the BAM file. Otherwise adapt here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Map transposon insertion sites
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%In the following loop: look for coordinates of each read and compares with
%%%%preceeding to figure out if it is the same transposon
ll=1;
for kk=1:17   %%%<-17 is number of chromosome (including mitochrondrial chromosome ChrMT)
    
    aa=BioMap(file,'SelectReference',infobam.SequenceDictionary(kk).SequenceName);
    start=aa.Start;
    flag=aa.Flag;
    seq=aa.Sequence;
    readlength=cellfun('length',seq); %Get the length of each nucleotide sequence
    clear seq
    
    %%%% correct the start site according to orientation and readlength(
    %%%% read length must be added to coordinate if read is in reverse
    %%%% orientation (i.e. flag==16))
    startdirect=start(flag==00);
    flagdirect=flag(flag==00);
    startindirect=start(flag==16)+uint32(readlength(flag==16));
    flagindirect=flag(flag==16);
    start2=[startdirect ; startindirect];
    flag2=[flagdirect ; flagindirect];
    clear flagdirect flagindirect startdirect startindirect
    
    %sort new coordinates
    [start2, sortmat]=sort(start2);
    flag2=flag2(sortmat);
    
    
    %process first read
    tncoordinates(ll,:)=[kk start2(1) flag2(1)]; %%
    tncoordinates=uint64(tncoordinates);
    mm=0; %counts how many reads there are with the same starting position and flag as the previous read.
    jj=1; %number of reads with unique starting positions or flags (i.e. reading orientations).
        for ii=2:size(start2,1) %loop starting from the second read (ii=2)
            %Save new information ONLY from the reads that do NOT have the
            %same starting position AND NOT the same flag.
            if abs(start2(ii)-start2(ii-1))<= 2 && flag2(ii)==flag2(ii-1)
                mm=mm+1;               
            else
                tncoordinates(ll,:)=[kk abs(mean([start2(ii-mm-1:ii-1)])) double(flag2(ii-1))];               
                mm=0;
                jj=jj+1;
                ll=ll+1; %==jj???
                              
            end
            readnumb(ll)=mm+1;
            
%             %This is just a counter that inform on progress as fraction of chromosome covered
%             if find([0:10000:10000000]==ii) %update only every 10000 steps
%                 infobam.SequenceDictionary(kk).SequenceName
%                 ii/size(start2,1)      
%             end
        end
            %%%%%%%LAST transposon on each chromosome is missed because
            %%%%%%%loop finishes too early to take it into account

    tnnumber(kk)=jj;
    infobam.SequenceDictionary(kk).SequenceName
    clear aa jj start start2 flag2 seq
end



%%%%%backup variables before transformation below
tncoordinates_copy=tncoordinates;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%prepare chromosomal features from gff file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%CDSs (stored as struct elements)
features.genes=find(strcmp(gff(:,3),'gene')); %get all genes from list
genes.chr=gff(features.genes,1); %get the chromosome where the genes belong to
genes.coordinates=cell2mat(gff(features.genes,[4 5])); %get the start and end basepair numbers of the genes

 genes.annotation=gff(find(strcmp(gff(:,3),'gene')),9) %% add annotation to genes. struct

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%put all features and  coordinate on one single huge concatenated chromosome
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tncoordinates_concat=tncoordinates_copy;
genes.coordinates_concat=genes.coordinates;
essential.coordinates_concat=essential.coordinates;

genes.chr_str=string(genes.chr); % to make genes.chr use same terminology as infobam.SequenceDictionary
for ic = 1:length(genes.chr_str)
if genes.chr_str(ic) == 'mt'
    genes.chr_str(ic) = 'Mito';
end
end

ll=0;
for ii=2:17
    ll=ll+infobam.SequenceDictionary(ii-1).SequenceLength; %Get the sequence length of the chromosome
    
    aa=find(tncoordinates_concat(:,1)==ii); %Find all indices of a specific chromosome.
    tncoordinates_concat(aa,2)=tncoordinates_concat(aa,2)+ll; %For all genes, the starting basepair is the basepair number on the chromosome plus the length of all previous chromosomes (concatenating the chromosomes).
    
    aa=find(strcmp(genes.chr_str,infobam.SequenceDictionary(ii).SequenceName)); %get all the genes from the current chromosome (represented by 'infobam.SequenceDictionary(ii).SequenceName').
    genes.coordinates_concat(aa,:)=genes.coordinates_concat(aa,:)+ll; %For both the start and end coordinates, add the length of the previous chromosome(s).
     
    aa=find(strcmp(essential.chr,infobam.SequenceDictionary(ii).SequenceName)); %get all the essential genes from the current chromosome.
    essential.coordinates_concat(aa,:)=essential.coordinates_concat(aa,:)+ll; %For both the start and end coordinates, add the length of the previous chromosome(s).
    
  
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%count number of transposon per gene minus end and beginning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


genelength = genes.coordinates(:,2)-genes.coordinates(:,1); %define gene length

for ii=1:length(genes.coordinates)
    %Find all reads (tncoordinates_copy) that start at a basepair number between the start and end of a gene (start_coor = gene.coordinates(ii,1) and end_coor = gene.coordinates(ii,2))
    xx=find(tncoordinates_concat(:,2)>=genes.coordinates_concat(ii,1)+genelength(ii)*0.1&tncoordinates_concat(:,2)<=genes.coordinates_concat(ii,2)-genelength(ii)*0.1);
    %determine how many reads there are per gene (index of each read is stored in xx). 
    tnpergeneminten(ii)=length(xx);
    readpergeneminten(ii)=sum(sum(readnumb(xx))-max(readnumb(xx)));
    rpgene_crude(ii)=sum(sum(readnumb(xx)));
end

%% Local read density regions of 100kb

wl = 100000; % define window length

roman = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "mt"];
u=0;
locallength = ones(length(genelength),1)*wl;

for k=1:length(chromosome_endpos)
    my_field = strcat('ch',num2str(k));     % create variable with all bp positions for every chromosome
    
    
    for j=u+1:u+sum(strcmp(roman(k),genes.chr)) % every gene of each chromosome
        if chromosome_endpos(k)<wl
            localreadnumb(j,1) = sum(readpos.(my_field)(1:chromosome_endpos(k)))/(sum(readpos.(my_field)(1:chromosome_endpos(k))~=0));
            locallength(j)=chromosome_endpos(17);
        elseif    genes.coordinates(j,1)-wl/2<1
            localreadnumb(j,1) = sum(readpos.(my_field)(1:wl))/(sum(readpos.(my_field)(1:wl)~=0));
        elseif genes.coordinates(j,1)+wl/2>chromosome_endpos(k)
            localreadnumb(j,1) = sum(readpos.(my_field)(chromosome_endpos(k)-wl:chromosome_endpos(k)))/(sum(readpos.(my_field)(chromosome_endpos(k)-wl:chromosome_endpos(k))~=0));
        else
            localreadnumb(j,1) = sum(readpos.(my_field)(genes.coordinates(j,1)-wl/2:genes.coordinates(j,1)+wl/2))/(sum(readpos.(my_field)(genes.coordinates(j,1)-wl/2:genes.coordinates(j,1)+wl/2)~=0));
        end
        u=j;
    end
    
end

localreaddensity=localreadnumb./locallength;
totalreaddensity=sum(readnumb)./sum(chromosome_endpos);

%% Calculates normalized read count per gene paper method

readpergene_norm_WTII_paper = readpergeneminten.*(totalreaddensity./localreaddensity').*(1./(genelength'.*sum(readnumb)));
readpergene_norm_WTII_paper(tnpergeneminten<3) = NaN;


%%
readratio_WTI_WT2_paper=readpergene_norm_WTI_paper./readpergene_norm_WTII_paper;
readratio_WTI_WT2_paper=readratio_WTI_WT2_paper';
%readratio_WTI_WT2(isnan(readratio_WTI_WT2)|isinf(readratio_WTI_WT2))=[];

%% plot

scatter(genelength,readratio_WTI_WT2_paper, '.')
set(gca,'yscale','log', 'Fontsize', 20)
xlabel('Genelength [bp]')
%xlim([0, 1100])
ylabel('Readratio WT1/WT2')
ylim([0.01, 100])

%% Calculates normalized read count per gene my method

%readpergene_norm_WTII = readpergeneminten_WT2./(tnpergeneminten_WT2.*sum(readnumb_WT2));
%readpergene_norm_WTII(tnpergeneminten_WT2<3) = NaN;

tnpergeneminten_WT1_WT2 = tnpergeneminten_WT1 + tnpergeneminten_WT2;

readpergene_norm_WT1 = (readpergeneminten_WT1./(tnpergeneminten_WT1-1)).*(length(tncoordinates_WT1)./sum(readnumb_WT1));
readpergene_norm_WT2 = (readpergeneminten_WT2./(tnpergeneminten_WT2-1)).*(length(tncoordinates_WT2)./sum(readnumb_WT2));
readpergene_norm_DplI = (readpergeneminten_DplI./(tnpergeneminten_DplI-1)).*(length(tncoordinates_DplI)./sum(readnumb_DplI));
readpergene_norm_DplIPsdII = (readpergeneminten_DplIPsdII./(tnpergeneminten_DplIPsdII-1)).*(length(tncoordinates_DplIPsdII)./sum(readnumb_DplIPsdII));
readpergene_norm_WT1_WT2 = ( readpergene_norm_WT1 + readpergene_norm_WT2 ) / 2;


readpergene_norm_WT1_WT2(tnpergeneminten_WT1<1) = NaN;
readpergene_norm_WT1_WT2(tnpergeneminten_WT2<1) = NaN;
readpergene_norm_DplI(tnpergeneminten_DplI<1) = NaN;
readpergene_norm_DplIPsdII(tnpergeneminten_DplIPsdII<1) = NaN;

%%

readratio_DplI_WT1WT2_mymethod=(readpergene_norm_DplI./readpergene_norm_WT1_WT2)';
readratio_DplIPsdII_WT1WT2_mymethod=(readpergene_norm_DplIPsdII./readpergene_norm_WT1_WT2)';
readratio_DplIPsdII_DplI_mymethod=(readpergene_norm_DplIPsdII./readpergene_norm_DplI)';

%%

readratio_WTI_WT2_mymethod_2=readpergene_norm_WTI_2./readpergene_norm_WTII_2;
readratio_WTI_WT2_mymethod_2=readratio_WTI_WT2_mymethod_2';

%% plot

scatter(tnpergene,readratio_WTI_WT2_mymethod_2, '.')
set(gca,'yscale','log', 'Fontsize', 20)
xlabel('Transposons per gene')
xlim([0, 500])
ylim([0.01,100])
ylabel('Readratio WT1/WT2')
