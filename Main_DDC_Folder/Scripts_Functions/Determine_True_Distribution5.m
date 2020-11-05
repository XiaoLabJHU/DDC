%This is a function that determines P_{blink} from the supporting material
%and defines the bins used within the likelihood calculation.

function [bins, D_Counts3, Total_No_Blink, Resolution]= Determine_Blinking_Distribution5(LocalizationsFinal, Frame_Information, Pre_A, Resolution)

%First we are going to go through and find the max of distance between locs
%in a image. We are also going to go through and find the min distance in
%an image
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


%Set the bins used in the likelihood calculation.
bins=[0:Resolution:D_maxf, Inf];

Total_Blink=[];
Total_No_Blink=[];

for i=1:length(LocalizationsFinal)
    
    disp(['Percent Done=',num2str(i/length(LocalizationsFinal)),' with STEP 1'])
    
    %Here we will gather the distance distrabutions from as many cells as
    %posible. This will allow us to define the joint probability
    %distribution more accuratly.
    
    Z2 = pdist([[1:length(Frame_Information{i})]'*0,Frame_Information{i}(:)]);
    
    D = (pdist(LocalizationsFinal{i}));
    
    D_Blink=D(Z2<Pre_A);
    
    D_No_Blink=D([(Z2>Pre_A).*(Pre_A*5>Z2)]==1);
    
    
    %Here we define the distrubions for each of the images, we then compair
    %them later.
    D_Counts=histcounts(D_Blink,bins,'Normalization','prob');
    
    Total_Blink=[Total_Blink; D_Counts];
    
    D_Counts2=histcounts(D_No_Blink,bins,'Normalization','prob');
    
    Total_No_Blink=[Total_No_Blink; D_Counts2];
    
end


%%
%Here we go through and scale the probablity distributions after 10 times
%the resolution of the experiment.

D_Counts=mean(Total_Blink,1);%histcounts(Total_Blink,bins,'Normalization','prob');

D_Counts2=mean(Total_No_Blink,1);%histcounts(Total_No_Blink,bins,'Normalization','prob');

%kk=find(D_Counts<D_Counts2,1);

D_Scale=sum(D_Counts(10:end))/sum(D_Counts2(10:end));

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
    disp('Warning: You may have defined the wrong resolution or you may have drifting, consider defining resolution yourself')
    D_Counts3(8:end)=0;
    D_Counts3=(D_Counts3)/sum(D_Counts3);
end

