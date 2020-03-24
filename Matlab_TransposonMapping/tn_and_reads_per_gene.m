%https://sites.google.com/site/satayusers/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Supplementary script 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[file, path]=uigetfile('*.bam');; %<-point to the BAM file
cd(path)

%%%load some variables
infobam=baminfo(file,'ScanDictionary',true);
load('yeastGFF.mat') %<- both yeastGFF.mat should be in the same folder as the BAM file. Otherwise adapt here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Map transposon insertion sites
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%Following loop: look for coordinates of each read and compares with
%%%%preceeding to figure out if it is the same transposon
ll=1
for kk=1:17   %%%<-17 is number of chromosome (including ChrMT)
    
    aa=BioMap(file,'SelectReference',infobam.SequenceDictionary(kk).SequenceName);
    start=aa.Start;
    flag=aa.Flag;
    seq=aa.Sequence;
    readlength=cellfun('length',seq); %Get the length of each nucleotide sequence
    clear seq
    
    %%%% correct the start site according to orientation and readlength(
    %%%% read length must be added to coordinate if read is in reverse
    %%%% orientation
    startdirect=start(flag==00);
    flagdirect=flag(flag==00);
    startindirect=start(flag==16)+uint32(readlength(flag==16));
    flagindirect=flag(flag==16);
    start2=[startdirect ; startindirect];
    flag2=[flagdirect ; flagindirect];
    clear flagdirect flagindirect startdirect startindirect
    
    %sort new coordinates
    [start2 sortmat]=sort(start2);
    flag2=flag2(sortmat);
    
    
    %process first read
    tncoordinates(ll,:)=[kk start2(1) flag2(1)]; %%
    tncoordinates=uint64(tncoordinates);
    mm=0; 
    jj=1;
        %loop
        for ii=2:size(start2,1)
            if abs(start2(ii)-start2(ii-1))<= 2 & flag2(ii)==flag2(ii-1);
                mm=mm+1;               
            else
                tncoordinates(ll,:)=[kk abs(mean([start2(ii-mm-1:ii-1)])) double(flag2(ii-1))];               
                mm=0;                
                jj=jj+1;
                ll=ll+1;
                              
            end
            readnumb(ll)=mm+1;
            
            %This is just a counter that inform on progress as fraction of chromosome covered
%             if find([0:10000:10000000]==ii)                     
%                 infobam.SequenceDictionary(kk).SequenceName     
%                 ii/size(start2,1)                                
%             end
        end
            %%%%%%%LAST transposon on each chromosome is missed because
            %%%%%%%loop finishes too early to take it into account

    tnnumber(kk)=jj
    clear aa jj start start2 flag2 seq
end



%%%%%backup variables before transformation below
tncoordinates_copy=tncoordinates;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%prepare chromosomal features from gff file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%CDSs
features.genes=find(strcmp(gff(:,3),'gene'));
genes.chr=gff(features.genes,1)
genes.coordinates=cell2mat(gff(features.genes,[4 5]));

%%%essential CDSs
features.essential=find(strcmp(gff(:,2),'YeastMine'));
essential.chr=gff(features.essential,1)
essential.coordinates=cell2mat(gff(features.essential,[4 5]));

%%%other features that can be optionally computed as well. Just uncomment
%%%what you are interested in.
% features.LTR_retrotransposon=find(strcmp(gff(:,3),'LTR_retrotransposon'));
% LTR_retrotransposon.chr=gff(features.LTR_retrotransposon,1);
% LTR_retrotransposon.coordinates=cell2mat(gff(features.LTR_retrotransposon,[4 5]));
% 
% features.tRNA=find(strcmp(gff(:,3),'tRNA'));
% tRNA.chr=gff(features.tRNA,1)
% tRNA.coordinates=cell2mat(gff(features.tRNA,[4 5]));
% 
% features.rRNA=find(strcmp(gff(:,3),'rRNA'));
% rRNA.chr=gff(features.rRNA,1)
% rRNA.coordinates=cell2mat(gff(features.rRNA,[4 5]));
% 
% features.centromere=find(strcmp(gff(:,3),'centromere'));
% centromere.chr=gff(features.centromere,1)
% centromere.coordinates=cell2mat(gff(features.centromere,[4 5]));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%put all features and  coordinate on one single huge concatenated chromosome
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ll=0
for ii=2:17
    ll=ll+infobam.SequenceDictionary(ii-1).SequenceLength;
    
    aa=find(tncoordinates_copy(:,1)==ii);
    tncoordinates_copy(aa,2)=tncoordinates_copy(aa,2)+ll;
    
    aa=find(strcmp(genes.chr,infobam.SequenceDictionary(ii).SequenceName));
    genes.coordinates(aa,:)=genes.coordinates(aa,:)+ll;
     
    aa=find(strcmp(essential.chr,infobam.SequenceDictionary(ii).SequenceName));
    essential.coordinates(aa,:)=essential.coordinates(aa,:)+ll;  
    
    %%%uncomment whatever is needed according to section above
%     aa=find(strcmp(promoters.chr,infobam.SequenceDictionary(ii).SequenceName));
%     promoters.coordinates(aa,:)=promoters.coordinates(aa,:)+ll;       
%        
%     aa=find(strcmp(LTR_retrotransposon.chr,infobam.SequenceDictionary(ii).SequenceName));
%     LTR_retrotransposon.coordinates(aa,:)=LTR_retrotransposon.coordinates(aa,:)+ll;
%     
%     aa=find(strcmp(tRNA.chr,infobam.SequenceDictionary(ii).SequenceName));
%     tRNA.coordinates(aa,:)=tRNA.coordinates(aa,:)+ll;
%     
%     aa=find(strcmp(rRNA.chr,infobam.SequenceDictionary(ii).SequenceName));
%     rRNA.coordinates(aa,:)=rRNA.coordinates(aa,:)+ll;
%     
%     aa=find(strcmp(centromere.chr,infobam.SequenceDictionary(ii).SequenceName));
%     centromere.coordinates(aa,:)=centromere.coordinates(aa,:)+ll;
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%count number of transposon per chromosomal features
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=1:length(genes.coordinates)
    xx=find(tncoordinates_copy(:,2)>=genes.coordinates(ii,1)&tncoordinates_copy(:,2)<=genes.coordinates(ii,2));
    tnpergene(ii)=length(xx);
    readpergene(ii)=sum(sum(readnumb(xx))-max(readnumb(xx)));
    rpgene_crude(ii)=sum(sum(readnumb(xx)));
end
for ii=1:length(essential.coordinates)
    xx=find(tncoordinates_copy(:,2)>=essential.coordinates(ii,1)&tncoordinates_copy(:,2)<=essential.coordinates(ii,2));
    tnperessential(ii)=length(xx);
    readperessential(ii)=sum(sum(readnumb(xx))-max(readnumb(xx)));
    rpessential_crude(ii)=sum(sum(readnumb(xx)));
end

%%%uncomment whatever is needed according to section above
% for ii=1:length(rRNA.coordinates)
%     xx=find(tncoordinates_copy(:,2)>=rRNA.coordinates(ii,1)&tncoordinates_copy(:,2)<=rRNA.coordinates(ii,2));
%     tnperrRNA(ii)=length(xx);
%     readperrRNA(ii)=sum(sum(readnumb(xx))-max(readnumb(xx)));
% end
% for ii=1:length(tRNA.coordinates)
%     xx=find(tncoordinates_copy(:,2)>=tRNA.coordinates(ii,1)&tncoordinates_copy(:,2)<=tRNA.coordinates(ii,2));
%     tnpertRNA(ii)=length(xx);
%     readpertRNA(ii)=sum(sum(readnumb(xx))-max(readnumb(xx)));
% end
% for ii=1:length(LTR_retrotransposon.coordinates)
%     xx=find(tncoordinates_copy(:,2)>=LTR_retrotransposon.coordinates(ii,1)&tncoordinates_copy(:,2)<=LTR_retrotransposon.coordinates(ii,2));
%     tnperLTR_retrotransposon(ii)=length(xx);
%     readperLTR_retrotransposon(ii)=sum(sum(readnumb(xx))-max(readnumb(xx)));
% end
% for ii=1:length(promoters.coordinates)
%     xx=find(tncoordinates_copy(:,2)>=promoters.coordinates(ii,1)&tncoordinates_copy(:,2)<=promoters.coordinates(ii,2));
%     tnperpromoter(ii)=length(xx);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Save .MAT file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename=[infobam.Filename,'.mat'];
save(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%generate and save bed file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%first, convert ChrMito into ChrM in infobam.SequenceDictionary(ii).SequenceName
infobam.SequenceDictionary(17).SequenceName='M';

filename=[infobam.Filename,'.bed']
fileID=fopen(filename,'wt');
str=['track name="' infobam.Filename '" useScore=1']
fprintf(fileID,str)
fprintf(fileID,'\n');

for ii=1:17
    Chr=['chr' char(infobam.SequenceDictionary(ii).SequenceName)];
    index=find(tncoordinates(:,1)==ii);
    for jj=1:length(index)
        str=sprintf('%1$s %2$s %3$s %4$s %5$s',Chr,num2str(tncoordinates(index(jj),2)),num2str(tncoordinates(index(jj),2)+1),'.',num2str(100+(readnumb(index(jj)))*20));
         fprintf(fileID,str);
         fprintf(fileID,'\n');
    end
ii
end
fclose(fileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%generate and save a table file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('names.mat')
filename=[infobam.Filename,'_pergene.txt']
fileID=fopen(filename,'wt');
str=['gene_name' char(9) 'number_of_transposon_per_gene' char(9) 'number_of_read_pergene']
fprintf(fileID,str)
fprintf(fileID,'\n');
for ii=1:length(names)
    str=sprintf('%s \t', char(names(ii)),num2str(tnpergene(ii)),num2str(readpergene(ii)));
    fprintf(fileID,str);
    fprintf(fileID,'\n');
end
fclose(fileID);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%write a wig file%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename=[infobam.Filename,'.wig']
fileID=fopen(filename,'wt');
str=['track type=wiggle_0 maxHeightPixels=60 name="' infobam.Filename '"']
fprintf(fileID,str)
fprintf(fileID,'\n');

%%%%find transposons mapping at same place, and sum them
%%%%

[~,cc]=unique(tncoordinates(:,[1 2]),'rows');
dd=setdiff(1:length(tncoordinates),cc);
tncoordinateswig=tncoordinates;
readnumbwig=readnumb;
readnumbwig(dd-1)=readnumb(dd-1)+readnumb(dd);
readnumbwig(dd)=[];
tncoordinateswig(dd,:)=[];
clear cc dd


for ii=1:16
    Chr=['chr' char(infobam.SequenceDictionary(ii).SequenceName)];
    str=['variableStep chrom=' Chr]
             fprintf(fileID,str);
         fprintf(fileID,'\n');
    index=find(tncoordinateswig(:,1)==ii);
    for jj=1:length(index)
        str=sprintf('%1$s %2$s %3$s %4$s %5$s',num2str(tncoordinateswig(index(jj),2)),num2str(readnumbwig(index(jj))));
         fprintf(fileID,str);
         fprintf(fileID,'\n');
    end
ii
end
fclose(fileID);
clear tncoordinateswig readnumbwig