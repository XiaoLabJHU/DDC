%This will be a function to determine where to set the bounds for the 'true
%distribution'. It will determine the maximum lifetime of the FP 'A' and
%determine how much data is availible at each frame difference, speeding up
%the analysis in a later section of the code. It will also gather data from
%every cell to determine "the average" blinking properties accross all
%cells.

function [bins, D_Counts3, Total_No_Blink]= Determine_True_Distribution5(LocalizationsFinal, Frame_Information, Pre_A, Resolution)

%This function will determine the maximum T_off for the particular PAF
%used and will determine what distances with which frame differences to
%include in the calculation of the true pairwise distance distribution.


%%
%First we are going to go through and find the max of distance
D_maxf=0;
for i=1:length(LocalizationsFinal)
    i/length(LocalizationsFinal)
    D = (pdist(LocalizationsFinal{i}));
    
    D_max=max(D);
    if D_max>D_maxf
        D_maxf=D_max;
    end
    
end


bins=[0:Resolution:D_maxf, Inf];



Total_Blink=[];
Total_No_Blink=[];
for i=1:length(LocalizationsFinal)
    
    %Here we need to get rid of the edge effect when calculating the
    %localizations pres
    %     IIs=find(LocalizationsFinal{i}(:,1)>(3*Resolution+min(LocalizationsFinal{i}(:,1))) & LocalizationsFinal{i}(:,1)<(-3*Resolution+max(LocalizationsFinal{i}(:,1))) & LocalizationsFinal{i}(:,2)>(3*Resolution+min(LocalizationsFinal{i}(:,2))) & LocalizationsFinal{i}(:,2)<(-3*Resolution+max(LocalizationsFinal{i}(:,2))));
    %
    %     Loc1=LocalizationsFinal{i}(IIs,:);
    %
    %
    %Here we will gather the distance distrabutions from as many cells as
    %posible. This will allow us to define the joint probability
    %distribution more accuratly.
    
    Z2 = pdist([[1:length(Frame_Information{i})]'*0,Frame_Information{i}(:)]);
    D = (pdist(LocalizationsFinal{i}));
    
    %     Dtest=pdist2(LocalizationsFinal{i},LocalizationsFinal{i});
    %     Dtest=Dtest(:);
    D_Blink=D(Z2<Pre_A);
    
    D_No_Blink=D([(Z2>Pre_A).*(Pre_A*5>Z2)]==1);
    
    
    D_Counts=histcounts(D_Blink,bins,'Normalization','prob');
    
    Total_Blink=[Total_Blink; D_Counts];
    
    D_Counts2=histcounts(D_No_Blink,bins,'Normalization','prob');
    
    Total_No_Blink=[Total_No_Blink; D_Counts2];
    
end


%%
%Instead of the old approach, we are going to use a kernal to estimate the
%probability distribution at Resolution/4 intervals.
%bins=[0:Resolution:max(Total_No_Blink)+Resolution*5, Inf];
D_Counts=mean( Total_Blink);%histcounts(Total_Blink,bins,'Normalization','prob');
%[D_Counts] = ksdensity(Total_Blink,bins,'Support','positive');
%[D_Counts2] = ksdensity(Total_No_Blink,bins,'Support','positive');
D_Counts2=mean(Total_No_Blink);%histcounts(Total_No_Blink,bins,'Normalization','prob');
% bins=[0:Resolution:max(Total_No_Blink)-Resolution, Inf];
%
%
% D_Counts=histcounts(Total_Blink,bins,'Normalization','prob');
% hold on
% D_Counts2=histcounts(Total_No_Blink,bins,'Normalization','prob');


kk=find(D_Counts<D_Counts2,1)
D_Scale=sum(D_Counts(10:end))/sum(D_Counts2(10:end));

D_Counts3=D_Counts-D_Counts2*D_Scale;
hold on 

plot(bins(1:end-1), D_Counts, bins(1:end-1), D_Counts2, bins(1:end-1), D_Counts2*(D_Scale))
drawnow
pause
%%
D_Counts3=(D_Counts3)/sum(D_Counts3);


%After 4 bins the probability to have a blink must be a monotonically
%decreasing function
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
%%
D_Counts3=(D_Counts3)/sum(D_Counts3);

%D_Counts3(ins)=D_Counts3(ins-1)/10;
D_Counts3=(D_Counts3)/sum(D_Counts3);
D_Counts=(D_Counts)/sum(D_Counts);
D_Counts2=(D_Counts2)/sum(D_Counts2);

plot(bins(1:end-1), D_Counts, bins(1:end-1), D_Counts2, bins(1:end-1), D_Counts3)

if sum(D_Counts3(8:end)>0)>1
    disp('Warning: You may have defined the wrong resolution or you may have drifting')
    D_Counts3(8:end)=0;
    D_Counts3=(D_Counts3)/sum(D_Counts3);
end
%%
clf
plot(bins(1:8), D_Counts3(1:8))
