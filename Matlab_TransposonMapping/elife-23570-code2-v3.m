%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear essentialdomain
for ii=1:length(genes.coordinates)	
    tnindex=find(tncoordinates_copy(:,2)>genes.coordinates(ii,1) & tncoordinates_copy(:,2)<genes.coordinates(ii,2));
    tns(ii).data=[genes.coordinates(ii,1); tncoordinates_copy(tnindex,2); genes.coordinates(ii,2)];
    tnnumber(ii)=length(tns(ii).data);
    genelength(ii)=genes.coordinates(ii,2)-genes.coordinates(ii,1);
end
 
mingap=200
maxlength=0.9
minlength=0.10
mintnnumber=20
skiptn=5 %minimum is 1 !
 
 
for ii=1:length(genes.coordinates)	
    if tnnumber(ii)<skiptn+1
    	tndif(ii)=0;
    else
        tndif(ii)=max(tns(ii).data(skiptn+1:length(tns(ii).data))-tns(ii).data(1:length(tns(ii).data)-skiptn));
    end
end
tndif=double(tndif);
essentialdomain=double(tndif).*tnnumber./genelength.^1.5;
essentialdomain(tndif<mingap | tndif./genelength>maxlength | tndif./genelength<minlength | tnnumber<mintnnumber)=0;
 
[dump index]=sort(essentialdomain,'descend');
clear dump
plot(essentialdomain(index),'.')
 
 
load('names.mat')
load('essentialgenesindex.mat')
 
 
names(index(1:200))
heatplot=zeros(length(essentialdomain),1);
heatplot(essentialgenes)=1;
figure(2)
for ii=1:6
    jj=ii*1000
    subplot(21,1,ii*3-2)
    image(rot90(256*heatplot(index(jj-999:jj))))
    colormap gray
end
figure(3)
for ii=1:6
    jj=ii*1000
    subplot(21,1,ii*3-1)
    image(50*(essentialdomain(index(jj-999:jj))))
    colormap jet
end
