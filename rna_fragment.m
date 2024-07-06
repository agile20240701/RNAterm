
% FUNCTION:  rna_fragment 

% Input parameters: 
%
% sequence: RNA sequence (from 5' to 3', a string containing C/U/A/G)
% charge:   precursor ion charge state (positive or negative + integer 1/2/3...)
% fiveprime: 5' terminal group ('OH' or 'PH', hydroxyl or phosphate group)
% threeprime: 3' terminal group ('OH' or 'OH', 'cP', cyclic phosphate included)
% modsite:  the site of modification in the sequence, 0: 5' terminus,1 to n
% modmass:  the monoisotopic mass of the modificaiton on the modsite
% xlswrt: 1,write all fragments into a .xls file named by the RNA sequence 
%         0,does not write fragments.

% Output parameters:
%
% fragments: the theoretical monoisotopic mass of all fragment ions
%            including multiply charged ions
% precursor_mz: the mass to charge ratio of the precursor ions.

% All Rights Reserved.
% By Dr. Ruixiang Sun, 
% National Insitute of Biological Sciences,Beijing,China
% Email: sunruixiang@nibs.ac.cn
% Built on: 2023-08-01
% Updated:  2024-07-04

function [fragments,precursor_mz]=rna_fragment(sequence,charge,fiveprime,threeprime,modsite,modmass,xlswrt)

%    for example
%    sequence='GACU'; 
%    charge=-4; 
%    fiveprime='OH';  or 'PH' 
%    threeprime='OH';  or 'PH','cP' 
%    modsite=-1; no modification
%    modmass=0.0; no modification
%    xlswrt=0; no xls file output

seqlen=length(sequence);  % seqlen: length of the input sequence.
sequence=upper(sequence); % uppercased to all letters input in sequence.
CUAG_comp=[length(find(sequence=='C')) length(find(sequence=='U')) length(find(sequence=='A')) length(find(sequence=='G'))];
%the number of C/U/A/G into CUAG_comp

% C H N O P  1C  2H            3N             4O            5P          
% monoisotopic mass of the elements of C H N O P into chnop_mass:
chnop_mass=[12.0 1.00782503214 14.00307400524 15.9949146221 30.973761998]';
%mass of a proton:
proton=1.00727645224;

%molecular formula of RNA four bases: C U A G, C1 H2 N3 O4
CUAG_bs_formula=[4 5 3 1; 4 4 2 2; 5 5 5 0; 5 5 5 1];  %atom number of C H N O
%nucleoside formula of C U A G
CUAG_ns_formula=[9 11 3 3; 9 10 2 4; 10 11 5 2; 10 11 5 3]; %atom number of C H N O for four nucleosides  

CUAG_ns_mass=CUAG_ns_formula*chnop_mass(1:4); %nucleoside mass
CUAG_bs_mass=CUAG_bs_formula*chnop_mass(1:4); %nucleobase mass

HO_mass=[0 1 0 1 0]*chnop_mass;    % mass of OH
H2O_mass=[0 2 0 1 0]*chnop_mass;   % mass of H2O
PO2H_mass=[0 -1 0 2 1]*chnop_mass; % mass of dehydrated metaphosphoric acid
HPO3_mass=[0 1 0 3 1]*chnop_mass;  % mass of metaphosphoric acid
HPO4_mass=[0 1 0 4 1]*chnop_mass;  % mass of HPO4
H3PO4_mass=[0 3 0 4 1]*chnop_mass; % mass of phosphoric acid

% mass of 5' terminal groups
if fiveprime=='OH'
    mass_5p=chnop_mass(2)+chnop_mass(4); %mass of 5' terminus for -OH
elseif fiveprime=='PH'
    mass_5p=chnop_mass(2)*2+chnop_mass(4)*4+chnop_mass(5); %mass of 5' terminus for phosphate group.
else
    error_num=5  %output error number '5'
end

% mass of 3' terminal groups
if threeprime=='OH'
    mass_3p=chnop_mass(2)+chnop_mass(4);
elseif threeprime=='PH'
    mass_3p=chnop_mass(2)*2+chnop_mass(4)*4+chnop_mass(5);
elseif threeprime=='cP' %cyclic phosphate group
    mass_3p=chnop_mass(4)*3+chnop_mass(5);
else
    error_num=3 %output error number '3'
end

%precursor mass, output precursor_mz
precursor=CUAG_comp*CUAG_ns_mass+(seqlen-1)*HPO4_mass+mass_5p+mass_3p+modmass; %five parts, mass summed up
if charge<0 %negatively charged case
    for i=1:-charge
        precursor_mz(i)=(precursor-i*proton)/i;
    end
else %positively charged case
    for i=1:charge
        precursor_mz(i)=(precursor+i*proton)/i;
    end
end

%mass of a1 and w1 ion (neutral mass)
first_loc=find(sequence(1)=='CUAG');
CUAG_comp_a=zeros(1,4);
CUAG_comp_a(first_loc)=1;
aion(1)=CUAG_ns_mass(first_loc)+mass_5p-chnop_mass(2); %mass of a1 ion
if modsite==0 | modsite==1
    aion(1)=aion(1)+modmass;  %a1 ion mass (here is neutral, proton mass will be added or substracted later)
end
%the number code to denote different types of RNA fragments:
% a/b/c/d/w/x/y/z: 11/12/13/14/15/16/17/18
seq_ions(1,1:13)=[11 1 CUAG_comp_a(1,1:4) aion(1) 0 0 0 0 first_loc 1]; %abcd/wxyz 11 12 13 14/15 16 17 18
seq_ions(1,7+first_loc)=aion(1)-CUAG_bs_mass(first_loc);

last_loc=find(sequence(seqlen)=='CUAG');
CUAG_comp_w=zeros(1,4);
CUAG_comp_w(last_loc)=1;
wion(1)=CUAG_ns_mass(last_loc)+mass_3p+chnop_mass(2)+HPO4_mass; %mass of w1 ion
if modsite==seqlen | modsite==seqlen+1
    wion(1)=wion(1)+modmass; % mass of modification
end
seq_ions(1+(seqlen-1)*4,1:13)=[15 seqlen CUAG_comp_w(1,1:4) wion(1) 0 0 0 0 last_loc seqlen];
seq_ions(1+(seqlen-1)*4,7+last_loc)=wion(1)-CUAG_bs_mass(last_loc);

%a1-base ion=aion(1)-CUAG_bs_mass(first_loc);

% mass of all a ions and w ions
for i=2:seqlen-1
    base_loc=find(sequence(i)=='CUAG');
    comp_ai=zeros(1,4);
    comp_ai(base_loc)=1;
    CUAG_comp_a(i,1:4)=CUAG_comp_a(i-1,1:4)+comp_ai;
    aion(i)=aion(i-1)+HPO4_mass+CUAG_ns_mass(base_loc);
    if modsite==i
        aion(i)=aion(i)+modmass;
    end
    seq_ions(i,1:13)=[11 1 CUAG_comp_a(i,1:4) aion(i) 0 0 0 0 base_loc i]; %!!! this line
    seq_ions(i,8:11)=aion(i)-CUAG_bs_mass'.*(CUAG_comp_a(i,1:4)>0);
    seq_ions(i,8:11)=seq_ions(i,8:11).*(CUAG_comp_a(i,1:4)>0);
    
    base_loc2=find(sequence(seqlen-i+1)=='CUAG');
    comp_wi=zeros(1,4);
    comp_wi(base_loc2)=1;
    CUAG_comp_w(i,1:4)=CUAG_comp_w(i-1,1:4)+comp_wi;
    wion(i)=wion(i-1)+HPO4_mass+CUAG_ns_mass(base_loc2);
    
    if modsite==seqlen-i+1
        wion(i)=wion(i)+modmass;
    end
    seq_ions(i+(seqlen-1)*4,1:13)=[15 seqlen-i+1 CUAG_comp_w(i,1:4) wion(i) 0 0 0 0 base_loc2 seqlen];
    seq_ions(i+(seqlen-1)*4,8:11)=wion(i)-CUAG_bs_mass'.*(CUAG_comp_w(i,1:4)>0);
    seq_ions(i+(seqlen-1)*4,8:11)=seq_ions(i+(seqlen-1)*4,8:11).*(CUAG_comp_w(i,1:4)>0);  

end

%b ions
seq_ions(seqlen:(seqlen-1)*2,1:13)=seq_ions(1:seqlen-1,1:13);
seq_ions(seqlen:(seqlen-1)*2,1)=12;
for i=seqlen:(seqlen-1)*2
    for j=7:11
        seq_ions(i,j)=seq_ions(i,j)+(seq_ions(i,j)>0)*H2O_mass;
    end
end

%c ions
seq_ions((seqlen-1)*2+1:(seqlen-1)*3,1:13)=seq_ions(1:seqlen-1,1:13);
seq_ions((seqlen-1)*2+1:(seqlen-1)*3,1)=13;
for i=(seqlen-1)*2+1:(seqlen-1)*3
    for j=7:11
        seq_ions(i,j)=seq_ions(i,j)+(seq_ions(i,j)>0)*HPO3_mass;
    end
end

%d ions
seq_ions((seqlen-1)*3+1:(seqlen-1)*4,1:13)=seq_ions(1:(seqlen-1),1:13);
seq_ions((seqlen-1)*3+1:(seqlen-1)*4,1)=14;
for i=(seqlen-1)*3+1:(seqlen-1)*4
    for j=7:11
        seq_ions(i,j)=seq_ions(i,j)+(seq_ions(i,j)>0)*H3PO4_mass;
    end
end

%x ions
seq_ions((seqlen-1)*5+1:(seqlen-1)*6,1:13)=seq_ions((seqlen-1)*4+1:(seqlen-1)*5,1:13);
seq_ions((seqlen-1)*5+1:(seqlen-1)*6,1)=16;
for i=(seqlen-1)*5+1:(seqlen-1)*6
    for j=7:11
        seq_ions(i,j)=seq_ions(i,j)-(seq_ions(i,j)>0)*H2O_mass;
    end
end

%y ions
seq_ions((seqlen-1)*6+1:(seqlen-1)*7,1:13)=seq_ions((seqlen-1)*4+1:(seqlen-1)*5,1:13);
seq_ions((seqlen-1)*6+1:(seqlen-1)*7,1)=17;
for i=(seqlen-1)*6+1:(seqlen-1)*7
    for j=7:11
        seq_ions(i,j)=seq_ions(i,j)-(seq_ions(i,j)>0)*HPO3_mass;
    end
end

%z ions
seq_ions((seqlen-1)*7+1:(seqlen-1)*8,1:13)=seq_ions((seqlen-1)*4+1:(seqlen-1)*5,1:13);
seq_ions((seqlen-1)*7+1:(seqlen-1)*8,1)=18;
for i=(seqlen-1)*7+1:(seqlen-1)*8
    for j=7:11
        seq_ions(i,j)=seq_ions(i,j)-(seq_ions(i,j)>0)*H3PO4_mass;
    end
end

%internal ions (two cleavages)
line=1;
internal=[0 H2O_mass 2*H2O_mass HPO3_mass H3PO4_mass H3PO4_mass+H2O_mass...
         2*HPO3_mass 2*HPO3_mass+H2O_mass 2*HPO3_mass+2*H2O_mass];
     % there are 9 masses for one fixed sequence.
for i=1:seqlen-2
    for j=1:seqlen-2-(i-1)
        i_ion(line,1)=j+1;
        i_ion(line,2)=j+i;
        inter_mass=0;
        for m=j+1:j+i
            inter_loc=find(sequence(m)=='CUAG');
            if m==j+1
                inter_mass=inter_mass+CUAG_ns_mass(inter_loc)-2*chnop_mass(2);
            else
                inter_mass=inter_mass+CUAG_ns_mass(inter_loc)+HPO4_mass;
            end
        end
        i_ion(line,3)=inter_mass;
        string=sequence(i_ion(line,1):i_ion(line,2));
        i_ion(line,5:8)=[length(find(string=='C')) length(find(string=='U')) length(find(string=='A')) length(find(string=='G'))];
        line=line+1;
    end
end

%number code of 31-39 for 9 different types of internal ions.
for i=1:line-1
    for j=1:9
        i_ions(j+(i-1)*9,1:8)=[i_ion(i,1:2) j+30 i_ion(i,3)+internal(j) i_ion(i,5:8)];
    end
end

%internal ions and their base loss
line_start=(seqlen-1)*8+1;
line_end=(seqlen-1)*8+size(i_ions,1);
seq_ions(line_start:line_end,1)=i_ions(:,3);
seq_ions(line_start:line_end,2)=i_ions(:,1);
seq_ions(line_start:line_end,3:6)=i_ions(:,5:8);
seq_ions(line_start:line_end,7)=i_ions(:,4);
seq_ions(line_start:line_end,12)=i_ions(:,1);
seq_ions(line_start:line_end,13)=i_ions(:,2);
for i=line_start:line_end
    seq_ions(i,8:11)=i_ions(i-line_start+1,4)-CUAG_bs_mass'.*(i_ions(i-line_start+1,5:8)>0);
    seq_ions(i,8:11)=seq_ions(i,8:11).*(i_ions(i-line_start+1,5:8)>0);
end

%Charged single base ions: C, U, A, G at low mass area type code=1
seq_ions(line_end+1,1:13)=[1 0 0 0 0 0 0 CUAG_bs_mass'.*(CUAG_comp(1:4)>0) 0 0];

%Precursor and its Base-loss: neutral mass,type code=2
seq_ions(line_end+2,1:2)=[2 1];
seq_ions(line_end+2,3:6)=CUAG_comp;
seq_ions(line_end+2,7)=precursor;
seq_ions(line_end+2,8:11)=precursor-CUAG_bs_mass'.*(CUAG_comp(1:4)>0);
seq_ions(line_end+2,8:11)=seq_ions(line_end+2,8:11).*(CUAG_comp(1:4)>0);
seq_ions(line_end+2,12:13)=[seqlen seqlen];

%Precursor-H20 and its Base-loss: neutral mass,type code=3
seq_ions(line_end+3,1:2)=[3 1];
seq_ions(line_end+3,3:6)=CUAG_comp;
seq_ions(line_end+3,7)=precursor-H2O_mass;
seq_ions(line_end+3,8:11)=precursor-CUAG_bs_mass'.*(CUAG_comp(1:4)>0);
seq_ions(line_end+3,8:11)=seq_ions(line_end+3,8:11).*(CUAG_comp(1:4)>0)-H2O_mass;
seq_ions(line_end+3,12:13)=[seqlen seqlen];

%Precursor-2*H20 and its Base-loss: neutral mass,type code=4
seq_ions(line_end+4,1:2)=[4 1];
seq_ions(line_end+4,3:6)=CUAG_comp;
seq_ions(line_end+4,7)=precursor-2*H2O_mass;
seq_ions(line_end+4,8:11)=precursor-CUAG_bs_mass'.*(CUAG_comp(1:4)>0);
seq_ions(line_end+4,8:11)=seq_ions(line_end+4,8:11).*(CUAG_comp(1:4)>0)-2*H2O_mass;
seq_ions(line_end+4,12:13)=[seqlen seqlen];

%%Precursor-CONH(43) and its Base-loss: neutral mass,type code=5
seq_ions(line_end+5,1:2)=[5 1];
seq_ions(line_end+5,3:6)=CUAG_comp;
seq_ions(line_end+5,7)=precursor-[1 1 1 1 0]*chnop_mass;
seq_ions(line_end+5,8:11)=precursor-CUAG_bs_mass'.*(CUAG_comp(1:4)>0);
seq_ions(line_end+5,8:11)=seq_ions(line_end+5,8:11).*(CUAG_comp(1:4)>0)-[1 1 1 1 0]*chnop_mass;
seq_ions(line_end+5,12:13)=[seqlen seqlen];

%Precursor-CONH(43)-H2O and its Base-loss: neutral mass,type code=6
seq_ions(line_end+6,1:2)=[6 1];
seq_ions(line_end+6,3:6)=CUAG_comp;
seq_ions(line_end+6,7)=precursor-[1 3 1 2 0]*chnop_mass;
seq_ions(line_end+6,8:11)=precursor-CUAG_bs_mass'.*(CUAG_comp(1:4)>0);
seq_ions(line_end+6,8:11)=seq_ions(line_end+6,8:11).*(CUAG_comp(1:4)>0)-[1 3 1 2 0]*chnop_mass;
seq_ions(line_end+6,12:13)=[seqlen seqlen];

% 158.925 (H2P2O6-Proton) type code=7 dehydrated pyrophosphate anion
seq_ions(line_end+7,1:13)=[7 0 0 0 0 0 159.93266 0 0 0 0 0 0];

if fiveprime=='PH' | threeprime=='PH'
%Precursor-97.9769 and its Base-loss: neutral mass,type code=8
seq_ions(line_end+8,1:2)=[8 1];
seq_ions(line_end+8,3:6)=CUAG_comp;
seq_ions(line_end+8,7)=precursor-[0 3 0 4 1]*chnop_mass;
seq_ions(line_end+8,8:11)=precursor-CUAG_bs_mass'.*(CUAG_comp(1:4)>0);
seq_ions(line_end+8,8:11)=seq_ions(line_end+8,8:11).*(CUAG_comp(1:4)>0)-[0 3 0 4 1]*chnop_mass;
seq_ions(line_end+8,12:13)=[seqlen seqlen];

%Precursor-97.9769-H2O and its Base-loss: neutral mass,type code=9
seq_ions(line_end+9,1:2)=[9 1];
seq_ions(line_end+9,3:6)=CUAG_comp;
seq_ions(line_end+9,7)=precursor-[0 5 0 5 1]*chnop_mass;
seq_ions(line_end+9,8:11)=precursor-CUAG_bs_mass'.*(CUAG_comp(1:4)>0);
seq_ions(line_end+9,8:11)=seq_ions(line_end+9,8:11).*(CUAG_comp(1:4)>0)-[0 5 0 5 1]*chnop_mass;
seq_ions(line_end+9,12:13)=[seqlen seqlen];

%Precursor-79.9663 and its Base-loss: neutral mass,type code=10
seq_ions(line_end+10,1:2)=[10 1];
seq_ions(line_end+10,3:6)=CUAG_comp;
seq_ions(line_end+10,7)=precursor-[0 1 0 3 1]*chnop_mass;
seq_ions(line_end+10,8:11)=precursor-CUAG_bs_mass'.*(CUAG_comp(1:4)>0);
seq_ions(line_end+10,8:11)=seq_ions(line_end+10,8:11).*(CUAG_comp(1:4)>0)-[0 1 0 3 1]*chnop_mass;
seq_ions(line_end+10,12:13)=[seqlen seqlen];

end

%all singly charged cation or anion mass
ionlen=size(seq_ions,1);
if charge>0
    for i=1:ionlen
        for j=7:11
            if seq_ions(i,j)>0
                seq_ions(i,j)=seq_ions(i,j)+proton;
            else
                seq_ions(i,j)=0;
            end
        end
    end
else
    for i=1:ionlen
        for j=7:11
            if seq_ions(i,j)>0
                seq_ions(i,j)=seq_ions(i,j)-proton;
            else
                seq_ions(i,j)=0;
            end
        end
    end
end

%sort all masses and then remove mass redundancy
[row col value]=find(seq_ions(:,7:11));
all_ion=[row col+6 value];
[ta tb]=sort(all_ion(:,3));
all_ions=all_ion(tb,:);

distinct_ions(1,1:6)=[all_ions(1,3) 1 seq_ions(all_ions(1,1),1) seq_ions(all_ions(1,1),2) ...
                      seq_ions(all_ions(1,1),13) all_ions(1,2)-7];

m=1;
n=7;
for i=2:size(all_ions,1)
    if abs(all_ions(i,3)-all_ions(i-1,3))<0.0001
        distinct_ions(m,n:n+3)=[seq_ions(all_ions(i,1),1) seq_ions(all_ions(i,1),2)...
                                seq_ions(all_ions(i,1),13) all_ions(i,2)-7];
        n=n+4;
        distinct_ions(m,2)=distinct_ions(m,2)+1;
    else
        n=7;
        m=m+1;
        distinct_ions(m,1:6)=[all_ions(i,3) 1 seq_ions(all_ions(i,1),1) seq_ions(all_ions(i,1),2) ...
                      seq_ions(all_ions(i,1),13) all_ions(i,2)-7];
    end
end

for i=1:size(distinct_ions,1)
   if distinct_ions(i,2)>1
      thisline=distinct_ions(i,3:distinct_ions(i,2)*4+2);
      [t1 t2]=sort(thisline(1:4:distinct_ions(i,2)*4-3));
      for j=1:distinct_ions(i,2)
          distinct_ions(i,(j-1)*4+3:(j-1)*4+6)=thisline((t2(j)-1)*4+1:(t2(j)-1)*4+4);
      end
   end
end

% for the same mass ions (isomers: code 11.18 for az; code 17.12 for yb; code 13.16 for cx; code 15.14 for wd
ion1118=find(distinct_ions(:,2)>1 & distinct_ions(:,6)==0 & distinct_ions(:,10)==0 & distinct_ions(:,3)==11 & distinct_ions(:,7)==18);
distinct_ions(ion1118,3)=11.18;

ion1217=find(distinct_ions(:,2)>1 & distinct_ions(:,6)==0 & distinct_ions(:,10)==0 & distinct_ions(:,3)==12 & distinct_ions(:,7)==17);
distinct_ions(ion1217,3)=17.12;
distinct_ions(ion1217,4:5)=distinct_ions(ion1217,8:9);

ion1316=find(distinct_ions(:,2)>1 & distinct_ions(:,6)==0 & distinct_ions(:,10)==0 & distinct_ions(:,3)==13 & distinct_ions(:,7)==16);
distinct_ions(ion1316,3)=13.16;

ion1415=find(distinct_ions(:,2)>1 & distinct_ions(:,6)==0 & distinct_ions(:,10)==0 & distinct_ions(:,3)==14 & distinct_ions(:,7)==15);
distinct_ions(ion1415,3)=15.14;
distinct_ions(ion1415,4:5)=distinct_ions(ion1415,8:9);

%multiple charged ion masses:
line_distinct=size(distinct_ions,1);
colum_distinct=size(distinct_ions,2);
charged_ions(:,1)=distinct_ions(:,1);
charged_ions(:,3:11)=distinct_ions(:,2:10);
charged_ions(:,2)=-1*ones(line_distinct,1);
if charge>0
   charged_ions(:,2)=ones(line_distinct,1);
end

if charge<-1 | charge>1
line_count=line_distinct+1;
for i=1:line_distinct
    base_number=charged_ions(i,6)-charged_ions(i,5)+1;
    if charged_ions(i,4)>1 & base_number>0 % three bases: maxcharge=2
        for j=2:min(base_number+1,abs(charge))  %max charge determination according to #nt and precursor charge
            if charge<-1
               charged_ions(line_count,1)=(charged_ions(i,1)-(j-1)*proton)/j;
               charged_ions(line_count,2)=-1*j;
            else
               charged_ions(line_count,1)=(charged_ions(i,1)+(j-1)*proton)/j;
               charged_ions(line_count,2)=j;  
            end
            charged_ions(line_count,3:11)=charged_ions(i,3:11);
            line_count=line_count+1;
        end
    end
end

[tmp1 tmp2]=sort(charged_ions(:,1));
charged_ions=charged_ions(tmp2,:);

end

count=1;
fragments(1,1:11)=charged_ions(1,1:11);
for i=2:size(charged_ions,1)
    if charged_ions(i,1)-charged_ions(i-1,1)>0.0002 & charged_ions(i,1)>100
        count=count+1;
        fragments(count,1:11)=charged_ions(i,1:11);
    elseif charged_ions(i,1)>100
        fragments(count,1:11)=charged_ions(i,1:11);
        count=count+1;
    end
end
if fragments(1,1)<100
   fragments=fragments(2:end,:); 
end


%write all fragments into the .xls file named by the RNA sequence 
if xlswrt==1
 xlswrite([sequence '_ions_SRX'],{'m/z','charge','#isomers','type_1','start_1','end_1','loss_1','type_2','start_2','end_2','loss_2'});
 xlswrite([sequence '_ions_SRX'],fragments,['A2:K' num2str(size(fragments,1)+1)]);
end

%#end of this function.