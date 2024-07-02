
% draw cip CID data of 8 raw files with different terminal groups
% built: 2023-11-08 by Sun Ruixiang.

clear all
%close all

load RNA8_CIP_CID

rnanum=2
charge=-1
fivet='OH';
threet='OH';

cidhcd=1; %1 cid; 2 hcd
nce_start=1;
nce_end=9;
Phos=1;  

tolerance=10; % ppm
intenthresh=0; % 2%above
mloss_thresh=1;
ionloss_thresh=1;
internalthresh=1;
xlswrt=0;

CUAG='CUAG';
rna40={'GUCA'
'GUCA'                  
'GUCA'                  
'GUCA'
'AUCGAUCG'                                  
'AUCGAUCG'                                                                  
'AUCGAUCG'                 
'AUCGAUCG'};

for i=1:8
   dist=rna_ions_0903(rna40{i},-1,fivet,threet,-1,0);
   precursormass(i)=dist(end,1);
end

for i=1:8
    negprecursor(1,i)=(precursormass(rnanum)+Phos*79.96633-(i-1)*1.00727645)/i;
    negprecursor(2,i)=negprecursor(1,i)+79.96633/i;
    negprecursor(3,i)=negprecursor(2,i)+79.96633/i;
end

% figure
% % set(gcf,'position',[10 130 1400 590])
% % set(gca,'position',[0.05 0.1 0.93 0.82])
% set(gcf,'position',[10 160 1450 570])  
% set(gca,'position',[0.05 0.1 0.93 0.82])
% specnum=find(RNA_header(:,10)==rnanum*2-cidhcd);
% stem(RNA_header(specnum,2),RNA_header(specnum,5),'linewidth',2,'markersize',10);
% h=line([negprecursor(1,1) negprecursor(1,1)],[0 max(RNA_header(specnum,5))]);
% set(h,'color','r','linewidth',3)
% h=line([negprecursor(2,1) negprecursor(2,1)],[0 max(RNA_header(specnum,5))]);
% set(h,'color','r','linewidth',2)
% h=line([negprecursor(3,1) negprecursor(3,1)],[0 max(RNA_header(specnum,5))]);
% set(h,'color','r','linewidth',1)
% 
% h=line([negprecursor(1,2) negprecursor(1,2)],[0 max(RNA_header(specnum,5))]);
% set(h,'color','c','linewidth',3)
% h=line([negprecursor(2,2) negprecursor(2,2)],[0 max(RNA_header(specnum,5))]);
% set(h,'color','c','linewidth',2)
% h=line([negprecursor(3,2) negprecursor(3,2)],[0 max(RNA_header(specnum,5))]);
% set(h,'color','c','linewidth',1)
% 
% h=line([negprecursor(1,3) negprecursor(1,3)],[0 max(RNA_header(specnum,5))]);
% set(h,'color','m','linewidth',3)
% h=line([negprecursor(2,3) negprecursor(2,3)],[0 max(RNA_header(specnum,5))]);
% set(h,'color','m','linewidth',2)
% h=line([negprecursor(3,3) negprecursor(3,3)],[0 max(RNA_header(specnum,5))]);
% set(h,'color','m','linewidth',1)
% 
% h=line([negprecursor(1,4) negprecursor(1,4)],[0 max(RNA_header(specnum,5))]);
% set(h,'color','b','linewidth',3)
% h=line([negprecursor(2,4) negprecursor(2,4)],[0 max(RNA_header(specnum,5))]);
% set(h,'color','b','linewidth',2)
% h=line([negprecursor(3,4) negprecursor(3,4)],[0 max(RNA_header(specnum,5))]);
% set(h,'color','b','linewidth',1)
% 
% h=line([negprecursor(1,5) negprecursor(1,5)],[0 max(RNA_header(specnum,5))/2]);
% set(h,'color','k','linewidth',3)
% h=line([negprecursor(2,5) negprecursor(2,5)],[0 max(RNA_header(specnum,5))/2]);
% set(h,'color','k','linewidth',2)
% h=line([negprecursor(3,5) negprecursor(3,5)],[0 max(RNA_header(specnum,5))/2]);
% set(h,'color','k','linewidth',1)
% 
% grid on


for mm=nce_start:nce_end
    NCE=mm*5+5;

scan=find(RNA_header(:,3)==charge & RNA_header(:,7)==NCE & RNA_header(:,10)==rnanum );
if isempty(scan)==1
    disp('No scan found!')
    break;
end
    
spec_scan=RNA_header(scan(1),1)

filename=filestr(rnanum).name;
msloc=find(filename=='-');
bsloc=find(filename=='_');
rna=filename(msloc(1)+1:msloc(2)-1);
rnaleng=length(rna);
fivep=filename(bsloc(1)+1:msloc(1)-1);
threep=filename(msloc(2)+1:bsloc(2)-1);

%compute the charged all fragment ions' mass/charge
charged_ions=rna_frag_240629(rna,charge,fivep,threep,-1,0,0);
%figure
%stem(charged_ions(:,1),90*ones(size(charged_ions,1)))

%match process
currentmz=RNA_mzinten(find(RNA_mzinten(:,3)==rnanum),1:2);
currentspec=currentmz(RNA_header(find(spec_scan==RNA_header(:,1) & RNA_header(:,10)==rnanum),8):RNA_header(find(spec_scan==RNA_header(:,1) & RNA_header(:,10)==rnanum),9),1:2);
minmz=min(currentspec(:,1));
relativeinten=currentspec(:,2)/max(currentspec(:,2))*100;
nce2=RNA_header(find(RNA_header(:,1)==spec_scan & RNA_header(:,10)==rnanum),7);
charge2=RNA_header(find(RNA_header(:,1)==spec_scan  & RNA_header(:,10)==rnanum),3);
if RNA_header(find(RNA_header(:,1)==spec_scan  & RNA_header(:,10)==rnanum),6)==1
    frag='CID';
else
    frag='HCD';
end
    
if rnanum<10
    figure('NumberTitle','off','Name',['#' filename(1:2)  '_' rna '_' frag '_' num2str(NCE)]);
else
    figure('NumberTitle','off','Name',['#' filename(1:2) '_' rna '_' frag '_' num2str(NCE)]);
end

h=stem(currentspec(:,1),relativeinten);
set(h,'Marker','none','linewidth',0.8,'color','k')
hold on

k=1;
for i=1:size(charged_ions,1)
    currentmass=charged_ions(i,1);
    ppm_value=(currentspec(:,1)-currentmass)/currentmass*10.^6;
    tp1=find(abs(ppm_value)<tolerance);
    if length(tp1)==1 & relativeinten(tp1)>intenthresh
        matched_peaks(k,1:16)=[tp1 currentspec(tp1,1) currentmass  ppm_value(tp1) currentspec(tp1,2) relativeinten(tp1) charged_ions(i,2:11)]; 
        k=k+1;
    elseif length(tp1)>1
        rnanum,i
        error='Multiple matched for the mass'
        mass=currentmass
    end
end

%figure
%tempabc=currentspec;
%stem(currentspec(find(relativeinten>intenthresh),1),relativeinten(find(relativeinten>intenthresh)))
%hold on
%stem(matched_peaks(find(matched_peaks(:,6)>intenthresh),2),matched_peaks(find(matched_peaks(:,6)>intenthresh),6),'r')

anno_rate(mm)=sum(matched_peaks(find(matched_peaks(:,6)>intenthresh),6))/sum(relativeinten(find(relativeinten>intenthresh)))*100;

if k==1
    disp('No matched peak!')
    break;
end

%found matched peaks
if size(matched_peaks,1)>0

    %find a-B ions and change their type 11 to 10.99
    for i=1:rnaleng-1
       base_loc(i)=find(upper(rna(i))=='CUAG');
       tp1=find(matched_peaks(:,9)==11 & matched_peaks(:,11)==i & matched_peaks(:,12)==base_loc(i));
       if length(tp1)>0
          matched_peaks(tp1,9)=10.99;
       end
    end
    
nterm=find(matched_peaks(:,12)==0 & matched_peaks(:,9)>10 & matched_peaks(:,9)<15);
abcd=matched_peaks(nterm,[2 4 6 7 9 11]);

for i=1:size(abcd,1)
    switch abcd(i,5) 
        case 11
          ions{i}=['a_{',num2str(abcd(i,6)) '}'];  
        case 12
          ions{i}=['b_{',num2str(abcd(i,6)) '}'];        
        case 13
          ions{i}=['c_{',num2str(abcd(i,6)) '}'];        
        case 14
          ions{i}=['d_{',num2str(abcd(i,6)) '}'];        
        case 11.18
          ions{i}=['a/z_{',num2str(abcd(i,6)) '}'];  
        case 13.16
          ions{i}=['c/x_{',num2str(abcd(i,6)) '}'];  
    end
  end
end

aions=find(abcd(:,5)==11);
if length(aions)>0
for i=1:size(aions,1)  
   if abcd(aions(i),4)<-1
     h1=text(abcd(aions(i),1)-10,abcd(aions(i),3)+3,ions(aions(i)));
     set(h1,'color','c','fontsize',10)
     h1=text(abcd(aions(i),1)+7,abcd(aions(i),3)+3,['^{' num2str(abcd(aions(i),4)) '}']);
     set(h1,'color','r','fontsize',10)
   else
     h1=text(abcd(aions(i),1)-4,abcd(aions(i),3)+1,ions(aions(i)));
     set(h1,'color','c','fontsize',10)
   end
end
h=stem(abcd(aions,1),abcd(aions,3));
set(h,'Marker','none','linewidth',1,'color','c','linewidth',2)
end

bions=find(abcd(:,5)==12);
if length(bions)>0
h1=text(abcd(bions,1)-4,abcd(bions,3)+1,ions(bions));
set(h1,'color',[0.01 0.66 0.62],'fontsize',12)
h=stem(abcd(bions,1),abcd(bions,3));
set(h,'Marker','none','linewidth',1,'color',[0.01 0.66 0.62],'linewidth',2)
end

cions=find(abcd(:,5)==13);
if length(cions)>0
for i=1:size(cions,1)
  if abcd(cions(i),4)<-1
    h1=text(abcd(cions(i),1)-5,abcd(cions(i),3)+1,ions(cions(i)));
    set(h1,'color',[0 0 1],'fontweight','bold','fontsize',15);
    h1=text(abcd(cions(i),1)+4,abcd(cions(i),3)+1.5,['^{' num2str(abcd(cions(i),4)) '}']);
    set(h1,'color','r','fontweight','bold','fontsize',12);
  else
    h1=text(abcd(cions(i),1)-5,abcd(cions(i),3)+1,ions(cions(i)));
    set(h1,'color',[0 0 1],'fontweight','bold','fontsize',15);
  end
end
h=stem(abcd(cions,1),abcd(cions,3));
set(h,'Marker','none','linewidth',1,'color',[0 0 1],'linewidth',3)
end

dions=find(abcd(:,5)==14);
if length(dions)>0
    for i=1:size(dions,1)
      if abcd(dions(i),4)<-1
         h1=text(abcd(dions(i),1)-5,abcd(dions(i),3)+1,ions(dions(i)));
         set(h1,'color',[0 0 1],'fontweight','normal','fontsize',10);
         h1=text(abcd(dions(i),1)+3,abcd(dions(i),3)+1.5,['^{' num2str(abcd(dions(i),4)) '}']);
         set(h1,'color','r','fontweight','bold','fontsize',8);
      else
         h1=text(abcd(dions,1)-5,abcd(dions,3)+1,ions(dions));
         set(h1,'color',[0.24 0.35 0.67],'fontsize',12)
      end
    end
h=stem(abcd(dions,1),abcd(dions,3));
set(h,'Marker','none','linewidth',1,'color',[0.24 0.35 0.67],'linewidth',2)
end

azions=find(abcd(:,5)==11.18);
if length(azions)>0
    for i=1:size(azions,1)
        if abcd(azions(i),4)<-1
            h1=text(abcd(azions(i),1)-10,abcd(azions(i),3)+1,ions(azions(i)));
            set(h1,'color','c','fontsize',10)            
            h1=text(abcd(azions(i),1)+8,abcd(azions(i),3)+1.5,[ '^{' num2str(abcd(azions(i),4)), '}']);
            set(h1,'color','c','fontsize',10)
        else
            h1=text(abcd(azions(i),1)-10,abcd(azions(i),3)+1,ions(azions(i)));
            set(h1,'color','c','fontsize',10)
        end
        h=stem(abcd(azions,1),abcd(azions,3));
        set(h,'Marker','none','color','c','linewidth',1)
    end
end

cxions=find(abcd(:,5)==13.16);
if length(cxions)>0
for i=1:size(cxions,1)
  if abcd(cxions(i),4)<-1
     h1=text(abcd(cxions(i),1)-10,abcd(cxions(i),3)+1,ions(cxions(i)));
     set(h1,'color','b','fontsize',10)
     h1=text(abcd(cxions(i),1)+12,abcd(cxions(i),3)+1.5,['^{' num2str(abcd(cxions(i),4)) '}']);
     set(h1,'color','b','fontsize',8)
  else
     h1=text(abcd(cxions,1)-10,abcd(cxions,3)+1,ions(cxions));
     set(h1,'color','b','fontsize',10)
  end
end
h=stem(abcd(cxions,1),abcd(cxions,3));
set(h,'Marker','none','color','b','linewidth',2)
end

cterm=find(matched_peaks(:,12)==0 & matched_peaks(:,9)>14 & matched_peaks(:,9)<19);
wxyz=matched_peaks(cterm,[2 4 6 7 9 10]);

for i=1:size(wxyz,1)
    switch wxyz(i,5) 
        case 15
          ions{i}=['w_{',num2str(rnaleng+1-wxyz(i,6)) '}'];  
        case 16
          ions{i}=['x_{',num2str(rnaleng+1-wxyz(i,6)) '}'];        
        case 17
          ions{i}=['y_{',num2str(rnaleng+1-wxyz(i,6)) '}'];         
        case 18
          ions{i}=['z_{',num2str(rnaleng+1-wxyz(i,6)) '}'];  
        case 17.12
          ions{i}=['y/b_{',num2str(rnaleng+1-wxyz(i,6)) '}'];            
        case 15.14
          ions{i}=['w/d_{',num2str(rnaleng+1-wxyz(i,6)) '}'];  
    end
end

wions=find(wxyz(:,5)==15);
if length(wions)>0
for i=1:size(wions,1)
  if wxyz(wions(i),4)<-1
    h1=text(wxyz(wions(i),1)-6,wxyz(wions(i),3)+1,ions(wions(i)));
    set(h1,'color',[0.69 0.09 0.12],'fontsize',12);
    h1=text(wxyz(wions(i),1)+6,wxyz(wions(i),3)+1.5,['^{' num2str(wxyz(wions(i),4)) '}']);
    set(h1,'color',[0 0 1],'fontweight','bold','fontsize',10);
  else
    h1=text(wxyz(wions(i),1)-6,wxyz(wions(i),3)+1,ions(wions(i)));
    set(h1,'color',[0.69 0.09 0.12],'fontsize',12);
  end
end
end

xions=find(wxyz(:,5)==16);
if length(xions)>0
h1=text(wxyz(xions,1)-5,wxyz(xions,3)+1,ions(xions));
set(h1,'color',[1 0 1],'fontsize',12)
end

yions=find(wxyz(:,5)==17);
if length(yions)>0
for i=1:size(yions,1)
  if wxyz(yions(i),4)<-1
    h1=text(wxyz(yions(i),1)-5,wxyz(yions(i),3)+2,ions(yions(i)));
    set(h1,'color',[1 0.38 0],'fontweight','bold','fontsize',15);
    h1=text(wxyz(yions(i),1)+8,wxyz(yions(i),3)+2,['^{' num2str(wxyz(yions(i),4)) '}']);
    set(h1,'color',[0 0 1],'fontweight','bold','fontsize',12);
  else
    h1=text(wxyz(yions(i),1)-5,wxyz(yions(i),3)+2,ions(yions(i)));
    set(h1,'color',[1 0.38 0],'fontweight','bold','fontsize',15);
  end
end
end

zions=find(wxyz(:,5)==18);
if length(zions)>0
h1=text(wxyz(zions,1)-5,wxyz(zions,3)+1,ions(zions));
set(h1,'color',[1 0.75 0.8],'fontsize',12)
end

wdions=find(wxyz(:,5)==15.14);
if length(wdions)>0
h1=text(wxyz(wdions,1)-10,wxyz(wdions,3)+1,ions(wdions));
set(h1,'color',[0.69 0.09 0.12],'fontsize',10)
end

ybions=find(wxyz(:,5)==17.12);
if length(ybions)>0
        for i=1:size(ybions,1)
        if wxyz(ybions(i),4)<-1
            h1=text(wxyz(ybions(i),1)-10,wxyz(ybions(i),3)+1,ions(ybions(i)));
            set(h1,'color',[1 0.38 0],'fontsize',10)
            h1=text(wxyz(ybions(i),1)+12,wxyz(ybions(i),3)+1.5,['^{' num2str(wxyz(ybions(i),4)) '}']); 
            set(h1,'color',[1 0.38 0],'fontsize',8)
        else
            h1=text(wxyz(ybions(i),1)-10,wxyz(ybions(i),3)+1,ions(ybions(i)));
            set(h1,'color',[1 0.38 0],'fontsize',10)
        end
    end
end

h=stem(wxyz(wions,1),wxyz(wions,3));
set(h,'Marker','none','linewidth',1,'color',[0.69 0.09 0.12],'linewidth',2)
h=stem(wxyz(wdions,1),wxyz(wdions,3));
set(h,'Marker','none','linewidth',1,'color',[0.69 0.09 0.12],'linewidth',2)

h=stem(wxyz(xions,1),wxyz(xions,3));
set(h,'Marker','none','linewidth',1,'color',[1 0 1],'linewidth',2)

h=stem(wxyz(yions,1),wxyz(yions,3));
set(h,'Marker','none','linewidth',1,'color',[1 0.38 0],'linewidth',3)
h=stem(wxyz(ybions,1),wxyz(ybions,3));
set(h,'Marker','none','linewidth',1,'color',[1 0.38 0],'linewidth',2)

h=stem(wxyz(zions,1),wxyz(zions,3));
set(h,'Marker','none','linewidth',1,'color',[1 0.75 0.80],'linewidth',2)

%internal and its base loss
internal=find(matched_peaks(:,9)>30 & matched_peaks(:,9)<40 & matched_peaks(:,6)>internalthresh);
internalions=matched_peaks(internal,[2 4 6 7 9 12]);

h=stem(internalions(:,1),internalions(:,3));
set(h,'Marker','none','linewidth',1,'color',[0 1 0.5],'linewidth',1.1)
for i=1:size(internalions,1)
  if internalions(i,4)<-1
     h1=text(internalions(i,1)-4,internalions(i,3)+2,['i^{' num2str(internalions(i,4)) '}']);
     set(h1,'color',[0 1 0.5],'fontweight','bold','fontsize',12)
  else
     h1=text(internalions(i,1)-2,internalions(i,3)+2,'i' );
     set(h1,'color',[0 1 0.5],'fontweight','bold','fontsize',12)
  end
end


loss=find(matched_peaks(:,9)<19 & matched_peaks(:,12)>0 & matched_peaks(:,6)>ionloss_thresh);
lossions=matched_peaks(loss,[2 4 6 7 9:12]);
h=stem(lossions(:,1),lossions(:,3));
set(h,'Marker','none','linewidth',1,'color',[0.6 0.2 0.98],'linewidth',1.1)
for i=1:size(lossions,1)
    if lossions(i,5)==1
            h1=text(lossions(i,1)-5,lossions(i,3)+2,[CUAG(lossions(i,8)) '^{-1}']);
            set(h1,'color',[0.6 0.2 0.98],'fontweight','bold','fontsize',10)                
    end
    
    if lossions(i,4)==-1
        if lossions(i,5)==10.99
              h1=text(lossions(i,1)-3,lossions(i,3)+1,['a_' num2str(lossions(i,7)) '-B(' CUAG(lossions(i,8)) ')']);
              set(h1,'color',[0.6 0.2 0.98],'fontweight','bold','fontsize',10)
        end
        if lossions(i,5)==11
              h1=text(lossions(i,1)-3,lossions(i,3)+1,['a_' num2str(lossions(i,7)) '-' CUAG(lossions(i,8)) ]);
              set(h1,'color',[0.6 0.2 0.98],'fontweight','normal','fontsize',8)
        end
        if lossions(i,5)==12
              h1=text(lossions(i,1)-3,lossions(i,3)+1,['b_' num2str(lossions(i,7)) '-' CUAG(lossions(i,8))]);
              set(h1,'color',[0.6 0.2 0.98],'fontweight','normal','fontsize',8)
        end
        if lossions(i,5)==13
              h1=text(lossions(i,1)-3,lossions(i,3)+1,['c_' num2str(lossions(i,7)) '-' CUAG(lossions(i,8))]);
              set(h1,'color',[0.6 0.2 0.98],'fontweight','normal','fontsize',8)
        end
        if lossions(i,5)==14
              h1=text(lossions(i,1)-3,lossions(i,3)+1,['d_' num2str(lossions(i,7)) '-' CUAG(lossions(i,8))]);
              set(h1,'color',[0.6 0.2 0.98],'fontweight','normal','fontsize',8)
        end
        
        if lossions(i,5)==15
              h1=text(lossions(i,1)-3,lossions(i,3)+1,['w_' num2str(lossions(i,7)-lossions(i,6)+1) '-' CUAG(lossions(i,8)) ]);
              set(h1,'color',[0.6 0.2 0.98],'fontweight','normal','fontsize',8)
        end
        if lossions(i,5)==16
              h1=text(lossions(i,1)-3,lossions(i,3)+1,['x_' num2str(lossions(i,7)-lossions(i,6)+1) '-' CUAG(lossions(i,8))]);
              set(h1,'color',[0.6 0.2 0.98],'fontweight','normal','fontsize',8)
        end
        if lossions(i,5)==17
              h1=text(lossions(i,1)-3,lossions(i,3)+1,['y_' num2str(lossions(i,7)-lossions(i,6)+1) '-' CUAG(lossions(i,8))]);
              set(h1,'color',[0.6 0.2 0.98],'fontweight','normal','fontsize',8)
        end
        if lossions(i,5)==18
              h1=text(lossions(i,1)-3,lossions(i,3)+1,['z_' num2str(lossions(i,7)-lossions(i,6)+1) '-' CUAG(lossions(i,8))]);
              set(h1,'color',[0.6 0.2 0.98],'fontweight','normal','fontsize',8)
        end
        
        if lossions(i,5)==2
              h1=text(lossions(i,1)-3,lossions(i,3)+1,['M-'  CUAG(lossions(i,8))]);
              set(h1,'color',[0.6 0.2 0.98],'fontweight','normal','fontsize',8)
        end
        if lossions(i,5)==3
              h1=text(lossions(i,1)-3,lossions(i,3)+1,['M-'  CUAG(lossions(i,8)) '-H2O']);
              set(h1,'color',[0.6 0.2 0.98],'fontweight','normal','fontsize',8)
        end
    end
      if lossions(i,4)<-1
          if lossions(i,5)==10.99
              h1=text(lossions(i,1)-3,lossions(i,3)+1,['[a_' num2str(lossions(i,7)) '-B(' CUAG(lossions(i,8)) ')]^{' num2str(lossions(i,4)) '}']);
              set(h1,'color',[0.6 0.2 0.98],'fontweight','bold','fontsize',10)
        end
        if lossions(i,5)==11
              h1=text(lossions(i,1)-3,lossions(i,3)+1,['[a_' num2str(lossions(i,7)) '-' CUAG(lossions(i,8)) ']^{' num2str(lossions(i,4)) '}']);
              set(h1,'color',[0.6 0.2 0.98],'fontweight','normal','fontsize',8)
        end
        if lossions(i,5)==12
              h1=text(lossions(i,1)-3,lossions(i,3)+1,['[b_' num2str(lossions(i,7)) '-' CUAG(lossions(i,8)) ']^{' num2str(lossions(i,4)) '}']);
              set(h1,'color',[0.6 0.2 0.98],'fontweight','normal','fontsize',8)
        end
        if lossions(i,5)==13
              h1=text(lossions(i,1)-3,lossions(i,3)+1,['[c_' num2str(lossions(i,7)) '-' CUAG(lossions(i,8)) ']^{' num2str(lossions(i,4)) '}']);
              set(h1,'color',[0.6 0.2 0.98],'fontweight','normal','fontsize',8)
        end
        if lossions(i,5)==14
              h1=text(lossions(i,1)-3,lossions(i,3)+1,['[d_' num2str(lossions(i,7)) '-' CUAG(lossions(i,8)) ']^{' num2str(lossions(i,4)) '}']);
              set(h1,'color',[0.6 0.2 0.98],'fontweight','normal','fontsize',8)
        end
        
        if lossions(i,5)==15
              h1=text(lossions(i,1)-3,lossions(i,3)+1,['[w_' num2str(lossions(i,7)-lossions(i,6)+1) '-' CUAG(lossions(i,8)) ']^{' num2str(lossions(i,4)) '}']);        
              set(h1,'color',[0.6 0.2 0.98],'fontweight','normal','fontsize',8)
        end
        if lossions(i,5)==16
              h1=text(lossions(i,1)-3,lossions(i,3)+1,['[x_' num2str(lossions(i,7)-lossions(i,6)+1) '-' CUAG(lossions(i,8)) ']^{' num2str(lossions(i,4)) '}']);
              set(h1,'color',[0.6 0.2 0.98],'fontweight','normal','fontsize',8)
        end
        if lossions(i,5)==17
              h1=text(lossions(i,1)-3,lossions(i,3)+1,['[y_' num2str(lossions(i,7)-lossions(i,6)+1) '-' CUAG(lossions(i,8)) ']^{' num2str(lossions(i,4)) '}']);
              set(h1,'color',[0.6 0.2 0.98],'fontweight','normal','fontsize',8)
        end
        if lossions(i,5)==18
              h1=text(lossions(i,1)-3,lossions(i,3)+1,['[z_' num2str(lossions(i,7)-lossions(i,6)+1) '-' CUAG(lossions(i,8)) ']^{' num2str(lossions(i,4)) '}']);
              set(h1,'color',[0.6 0.2 0.98],'fontweight','normal','fontsize',8)
        end
        
        if lossions(i,5)==2
              h1=text(lossions(i,1)-3,lossions(i,3)+1,['[M-'  CUAG(lossions(i,8)) ']^{' num2str(lossions(i,4)) '}']);
              set(h1,'color',[0.6 0.2 0.98],'fontweight','normal','fontsize',8)
        end
        if lossions(i,5)==3
              h1=text(lossions(i,1)-3,lossions(i,3)+1,['[M-'  CUAG(lossions(i,8)) '-H2O]^{' num2str(lossions(i,4)) '}']);
              set(h1,'color',[0.6 0.2 0.98],'fontweight','normal','fontsize',8)
        end
        if lossions(i,5)==5
              h1=text(lossions(i,1)-3,lossions(i,3)+1,['[M-'  CUAG(lossions(i,8)) '-CONH-H2O]^{' num2str(lossions(i,4)) '}']);
              set(h1,'color',[0.6 0.2 0.98],'fontweight','normal','fontsize',8)
        end
      end
 
 end

%a-B draw



precursor=find(matched_peaks(:,9)<10 & matched_peaks(:,12)==0 & matched_peaks(:,6)>intenthresh);
preions=matched_peaks(precursor,[2 4 6 7 9 11]);
h=stem(preions(:,1),preions(:,3));
set(h,'Marker','none','linewidth',1,'color',[0 0 0],'linewidth',1.1)
for i=1:size(preions,1)
    if preions(i,4)<-1 & preions(i,3)>intenthresh
        switch preions(i,5)
            case 2
                h1=text(preions(i,1)-5,preions(i,3)+2,['M^{' num2str(preions(i,4)) '}']);
                set(h1,'color',[0 0 0],'fontweight','normal','fontsize',8)
            case 3
                h1=text(preions(i,1)-5,preions(i,3)+2,['[M-H2O]^{' num2str(preions(i,4)) '}']);
                set(h1,'color',[0 0 0],'fontweight','normal','fontsize',8)
            case 5  
                h1=text(preions(i,1)-5,preions(i,3)+2,['[M-CONH]^{' num2str(preions(i,4)) '}']);
                set(h1,'color',[0 0 0],'fontweight','normal','fontsize',8)
            case 8
                h1=text(preions(i,1)-5,preions(i,3)+2,['[M-H3PO4]^{' num2str(preions(i,4)) '}']);
                set(h1,'color',[0 0 0],'fontweight','normal','fontsize',8)
        end
    elseif preions(i,3)>intenthresh
        switch preions(i,5)
            case 2
                h1=text(preions(i,1)-5,preions(i,3)+1.5,'M');
                set(h1,'color',[0 0 0],'fontweight','normal','fontsize',8)
            case 3
                h1=text(preions(i,1)-5,preions(i,3)+1.5,'M-H2O');
                set(h1,'color',[0 0 0],'fontweight','normal','fontsize',8)                
            case 5
                h1=text(preions(i,1)-5,preions(i,3)+1.5,'M-CONH');
                set(h1,'color',[0 0 0],'fontweight','normal','fontsize',8)
            case 8
                h1=text(preions(i,1)-5,preions(i,3)+1.5,'M-H3PO4');
                set(h1,'color',[0 0 0],'fontweight','normal','fontsize',8)
        end
    end
end

line_matched=size(matched_peaks,1);
%178.235 peak
currentmass=178.235;
ppm_value=(currentspec(:,1)-currentmass)/currentmass*10.^6;
tp1=find(abs(ppm_value)<10);
if length(tp1)==1
    h=stem(currentspec(tp1(1),1),max(relativeinten(tp1)));
    set(h,'Marker','none','linewidth',0.5,'color',[0 0 0],'linewidth',0.5)
    h=text(currentspec(tp1(1),1)-5,max(relativeinten(tp1))+1.5,'CM');
    set(h,'color',[0 0 0],'fontweight','bold','fontsize',6)
    matched_peaks(line_matched+1,1:9)=[tp1 currentspec(tp1,1) 178.235 ppm_value(tp1) currentspec(tp1,2) relativeinten(tp1) -1 1 7];
else if length(tp1)>1
        error='Multiple matched for the mass 178.235'
        mass=currentspec(tp1,1:2)
    end
end

%158.9254 peak
line_matched=size(matched_peaks,1);
currentmass=158.9254;
ppm_value=(currentspec(:,1)-currentmass)/currentmass*10.^6;
tp1=find(abs(ppm_value)<10);
if length(tp1)==1
    h=stem(currentspec(tp1(1),1),max(relativeinten(tp1)));
    max(relativeinten(tp1))
    set(h,'Marker','none','linewidth',3,'color', 'r','linewidth',3)
    h=text(currentspec(tp1(1),1)-4,max(relativeinten(tp1))+1,'PP');
    set(h,'color','r','fontweight','bold','fontsize',14)
    matched_peaks(line_matched+1,1:9)=[tp1 currentspec(tp1(1),1) 158.9254 ppm_value(tp1(1)) currentspec(tp1,2) relativeinten(tp1) -1 1 7];
else if length(tp1)>1
        error='Multiple matched for the mass 158.9254'
        mass=currentspec(tp1,1:2)
    end
end

%98 and 98+H20 loss peaks
precursor_98=find(matched_peaks(:,9)>7 & matched_peaks(:,9)<10  & matched_peaks(:,6)>internalthresh);
preions_98=matched_peaks(precursor_98,[2 4 6 7 9 11]);
h=stem(preions_98(:,1),preions_98(:,3));
set(h,'Marker','none','linewidth',1,'color',[0 0 0],'linewidth',1.1)
for i=1:size(preions_98,1)
    if preions_98(i,4)<-1
        h1=text(preions_98(i,1)-5,preions_98(i,3)+2,['-P^{' num2str(preions_98(i,4)) '}']);
        set(h1,'color',[0 0 0],'fontweight','bold','fontsize',10)
    else
        h1=text(preions_98(i,1)-5,preions_98(i,3)+1.5,'-P');
        set(h1,'color',[0 0 0],'fontweight','bold','fontsize',10)
    end
end

xlabel('\it{m/z}')
ylabel('Relative Abundance (%)')
set(gca,'fontsize',12,'fontweight','bold')
%text(min(currentspec(:,1))-48,102,num2str(max(currentspec(:,2)),'%1.1e'),'fontsize',8)
text(5,102,num2str(max(currentspec(:,2)),'%1.1e'),'fontsize',8)
%grid on

% title([filename(1:2) '-' rna '(' num2str(rnaleng) ')-' frag ' NCE=' num2str(nce2) ' scan=' num2str(spec_scan) ' charge=' num2str(charge2) ' m/z=' num2str(charged_ions(end,1)) ' 5p=' fivep ' 3p=' threep])
% text(min(currentspec(:,1)),95,rna,'fontweight','bold','fontsize',20)
set(gcf,'position',[10 160 1450 570])  
set(gca,'position',[0.05 0.1 0.93 0.82])
axis([min(currentspec(:,1))-50 max(currentspec(:,1))+50 0 105])
%axis([100 700 0 105])

last=currentspec;
distlast=charged_ions;
clear currentmass i k ppm_value tp1 currentmz currentspec tp1 nterm abcd
clear aions bions cions dions wions xions yions zions 

if nce_end-nce_start>0
  clear matched_peaks
end

end

ax=gca;
ax.Box='off';
%set(gca,'xlim',[90 1220],'xtick',100:100:1200,'ytick',0:20:100);

if fivep=='PH'
    fivep='P';
end
if threep=='PH'
    threep='P';
end

%text(minmz,98,[fivep '-' rna '-' threep],'fontsize',15,'fontweight','bold');
%text(minmz,94,['\itz' '=' num2str(-charge) '^{-} HCD-50'],'fontsize',15,'fontweight','bold');
text(105,98,[fivep '-' rna '-' threep],'fontsize',15,'fontweight','bold');
text(102,94,['\itz' '=' num2str(-charge) '^{-} HCD-' num2str(NCE)],'fontsize',15,'fontweight','bold');


if xlswrt==1
 xlswrite([num2str(rnanum) '_' rna '_' frag num2str(charge) '_nce' num2str(NCE) '_scan_' num2str(spec_scan)],{'pos.','exp_mass','theo_mass','ppm','Inten','RelaInten','charge','#repeats','type_1','start_1','end_1','loss_1','type_2','start_2','end_2','loss_2'});
 xlswrite([num2str(rnanum) '_' rna '_' frag num2str(charge) '_nce' num2str(NCE) '_scan_' num2str(spec_scan)],matched_peaks,['A2:P' num2str(size(matched_peaks,1)+1)]);
end

%matched_peaks(find(matched_peaks(:,9)==7),6);
%PP_ex_th_ppm_RA=matched_peaks(end,2:6)%find(matched_peaks(:,2)<158.926 & matched_peaks(:,2)>158.925),6)

%text(162,98,['Exp= ' num2str(PP_ex_th_ppm_RA(1))],'fontsize',12,'fontweight','bold','color','r');
%text(162,95,['Theo= ' num2str(PP_ex_th_ppm_RA(2))],'fontsize',12,'fontweight','bold','color','r');
%text(162,92,['RA= ' num2str(PP_ex_th_ppm_RA(5)) '%'],'fontsize',12,'fontweight','bold','color','r');

anno_rate
figure
plot(10:5:50,anno_rate,'r+-')