%This is going to be a new script that will determine the nessary
%resolution of the bins for DDC. Making DDC more user friendly and
%availible to the public. 
clear
%First we are going to load in the data.

Location=pwd;
addpath([Location, '/Scripts_Functions'])

% User input for data file
[filename, pathname] = uigetfile({'*.mat'}, 'Select .mat file');
if ( filename == 0 )
    disp('Error! No (or wrong) file selected!')
    return
end

% Load all localizations from data structure file
full_filename = [ pathname, filename ];
load(full_filename,'Frame_Information','filename','TrueLocalizations','LocalizationsFinal');

%Here is where you are going to have to put in your determine frame where
%the steady state is reached, aka N!!!!!!!!!!!!!!!!!!!!
N_f=200;


%%
%First we are going to go through and find the max of distance between locs
%in a image. We are also going to go through and find the min distance in
%an image
Pre_A=N_f;
D_maxf=0;
D_maxm=Inf;
for i=1:length(LocalizationsFinal)
    i/length(LocalizationsFinal);
    D = (pdist(LocalizationsFinal{i}));
    
    D_max=max(D);
    if D_max>D_maxf
        D_maxf=D_max;
    end
    
    if D_max<D_maxm
        D_maxm=D_max;
    end
    
end


%%



%Here we are going to do an iterative process to determine what the proper
%bin size should be. We will start with a very large resolution for the
%image and then work our way down in 10nm segments. The smallest we are
%going to go is going to be 10nm. Once the value in the first bin is not
%the largest, that is when DDC assumes that the bin size is correct.

%Note: you will still want to look at the distribution to make sure that 

disp(['Looking for the correct bin resolution, this should only take a few min']) 
    
for Resolution=150:-10:10
disp(['Still working making prog, checking res of ',num2str(Resolution)])
%Set the bins used in the likelihood calculation.
bins=[0:Resolution:D_maxf, Inf];

Total_Blink=[];

Total_No_Blink=[];

for i=1:length(LocalizationsFinal)
    
    
    Z2 = pdist([[1:length(Frame_Information{i})]'*0,Frame_Information{i}(:)]);
    
    D = (pdist(LocalizationsFinal{i}));
    
    D_Blink=D(Z2<Pre_A);
    
    D_No_Blink=D([(Z2>Pre_A)]==1);
    
    %Here we define the distrubions for each of the images, we then compair
    %them later.
    D_Counts=histcounts(D_Blink,bins,'Normalization','prob');
    
    Total_Blink=[Total_Blink; D_Counts];
    
    D_Counts2=histcounts(D_No_Blink,bins,'Normalization','prob');
    
    Total_No_Blink=[Total_No_Blink; D_Counts2];
    
end

D_Counts=mean(Total_Blink,1);%histcounts(Total_Blink,bins,'Normalization','prob');

D_Counts2=mean(Total_No_Blink,1);%histcounts(Total_No_Blink,bins,'Normalization','prob');

for lsne=1:10000
    if lsne*Resolution>1000
        break
    end
end

D_Scale=sum(D_Counts(lsne:end))/sum(D_Counts2(lsne:end));

D_Counts3=D_Counts-D_Counts2*D_Scale;

D_Counts3=(D_Counts3)/sum(D_Counts3);

%Finally we clean up the distributions a little bit so that the logic is
%consistent with how the distribution should behave
good=1;
ins=0;
for i=4:length(D_Counts3)-1
    if D_Counts3(i)>0 && D_Counts3(i+1)<D_Counts3(i) && good==1
        
    else
        if ins==0
            ins=i;
        end
        good=0;
        D_Counts3(i)=0;
    end
end

D_Counts3(D_Counts3<0)=0;
D_Counts3=(D_Counts3)/sum(D_Counts3);
D_Counts3=(D_Counts3)/sum(D_Counts3);

%Anything greater than 8 times the resolution of the experiment is
%considered noise and will be eliminated from the probability distribution.
%We then warn the user if this is the case as they may want to consider a
%lower resolution.

if sum(D_Counts3(8:end)>0)>1
    disp('Warning: Eliminating Noise for higher bins')
    D_Counts3(8:end)=0;
    D_Counts3=(D_Counts3)/sum(D_Counts3);
end

% figure(1)
% plot(D_Counts3,'DisplayName','Total_No_Blink')
[I1,I2]=max(D_Counts3);


%Here is the criterion where the bin resolution is determined!! 
if I2~=1 

    Resolution=Resolution+10;
    break
end

end

disp(['The Resolution you should use for DDC is=',num2str(Resolution)])