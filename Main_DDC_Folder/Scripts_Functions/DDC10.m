
%This is the main section of the algorithm that will eliminate blinking in
%PALM or STORM data. Please see the user guide to understand the logic of
%each section.

%Christopher Herrick Bohrer 2018, Lab of Jie Xiao and Elijah Roberts

function [Final_Localizations_Blinking_Corrected, Final_Frame_Blinking_Corrected, LikHood, Score,  Dscale_store]=DDC10(ksu, LocalizationsFinal, Frame_Information, A, Resolution, TrueLocalizations, bins, Distribution_for_Blink, Total_No_Blink)


%First we will organize the localizations based off of the frame, so that
%they are in order.
[B, Inds]=sort(Frame_Information);
Frame_Information=Frame_Information(Inds);
LocalizationsFinal=LocalizationsFinal(Inds,:);

DistanceMatrix={};
for ijk=1:max(Frame_Information)+1
    DistanceMatrix{ijk}=[];
end

[DistanceMatrix]=Make_Distance_Matrix_Individual(LocalizationsFinal, DistanceMatrix, Frame_Information, max(Frame_Information));
DistanceControl2=[];

Proportion_In_Blinking_Range=[];
%Here we generate the True distribution using the distances that
for i=1:max(Frame_Information)
    
    if length(DistanceMatrix{i}(:))>100 && i>A
        DistanceControl2=[DistanceControl2; DistanceMatrix{i}(:)];
    end
    
    Proportion_In_Blinking_Range(end+1)=sum(DistanceMatrix{i}(:)<200)/length(DistanceMatrix{i}(:));
    
end

if length(DistanceControl2)<2000
    disp('WARNING: You do not have enough data to define a true distribution, if you have split an image consider combining multiple images!!')
end

True_Distribuiton=((histcounts(DistanceControl2,bins,'Normalization','prob')));
%[True_Distribuiton2] = ksdensity(DistanceControl2,bins,'Support','positive');

%Lets just convert it into probabilities
True_Distribuiton=True_Distribuiton/sum(True_Distribuiton);

if length(Total_No_Blink)>1
    D_Overall = (pdist(LocalizationsFinal));
    True_Distribuiton=((histcounts(Total_No_Blink(Total_No_Blink<max(D_Overall )),bins,'Normalization','prob')));
end
%Here we will determine a rough approximation that a localization is a
%blink, we will also generate the probability distributions that a distance
%belonging to a blink will have, given a frame separation an certain
%distance.

[Deviation_in_Probability, Prob_Distributions, Dscale_store]=Determine_Deviation_in_Probability8(True_Distribuiton,Distribution_for_Blink, bins, A, LocalizationsFinal, Frame_Information, [], []);

Prob_Distributions=log(Prob_Distributions);

%
% Deviation_in_Probability(isnan(Deviation_in_Probability(:,:)))=0;
% Deviation_in_Probability(isinf(Deviation_in_Probability(:,:)))=0;
%
% CalculationData=Deviation_in_Probability;
% CalculationData=CalculationData.*(CalculationData>0);
% CalculationData(CalculationData>1)=0;
% CalculationData(isnan(CalculationData(:,:)))=0;
% %
% Prob_of_Blink_Calculated2=Determine_Prob_of_Blinking2(LocalizationsFinal, A, Frame_Information, CalculationData, Resolution);
%
% loc=LocalizationsFinal(Prob_of_Blink_Calculated2==0,:);
% fram=Frame_Information(Prob_of_Blink_Calculated2==0);
%
% [~, Prob_Distributions]=Determine_Deviation_in_Probability78(True_Distribuiton,Distribution_for_Blink, bins, A, LocalizationsFinal, Frame_Information, loc, fram);
% Prob_Distributions=log(Prob_Distributions);

step_total=1000;

%True localizations if availible
if length(TrueLocalizations)>1
    True=[];
    True_frame=[];
    New_Loc=[];
    New_Frame=[];
    vvseo=zeros(1,length(LocalizationsFinal(:,1)));
    for kksd=1:length(LocalizationsFinal(:,1))
        
        if (LocalizationsFinal(kksd,1)<-1000) && sum(LocalizationsFinal(kksd,1)==TrueLocalizations(:,1))>0
            New_Frame=[New_Frame; Frame_Information(kksd)];
            New_Loc=[New_Loc; LocalizationsFinal(kksd,:)];
            
        elseif (LocalizationsFinal(kksd,1)<-1500)
            New_Frame=[New_Frame; Frame_Information(kksd)];
            New_Loc=[New_Loc; LocalizationsFinal(kksd,:)];
        end
        
        if sum(find(LocalizationsFinal(kksd,1)==TrueLocalizations(:,1)))>0
            True=[True; LocalizationsFinal(kksd,:)];
            True_frame=[True_frame; Frame_Information(kksd)];
            vvseo(kksd)=1;
        else
            
        end
    end
    loc1=True;
    loc=True;
    [denstiy_map1]=PALMplot(loc, loc1);
end
%
%  Prob_of_Blink_Calculated2=Determine_Prob_of_Blinking2_only_true(LocalizationsFinal, A, Frame_Information, CalculationData, Resolution, True);
%
%
% % %here we will calculate the amount of background prob will be generated
% % %just based off of the density and frame of each molecule
% Adjust_for_Background2
% %
% % %Prob of being a blink and not having it come from the random positioning
% % %of the molecule
%  Prob_of_Blink_Calculated2=Prob_of_Blink_Calculated2.*(1-final_adjust_prob);

%valo=Prob_of_Blink_Calculated2./final_adjust_prob;
%valo(valo>50)=50;
%We will maximize Final_lik
Final_lik=-Inf;
lik=999999;

%Calculate the pairwise distances and the overall frame differences
Z2 = pdist([[1:length(Frame_Information)]'*0,Frame_Information(:)]);

% for a stupid reason some people round their localizations. 
LocalizationsFinaltt=LocalizationsFinal;
good=0;
while good==0
    
LocalizationsFinaltt=LocalizationsFinal+rand(size(LocalizationsFinal));
D_Overall = (pdist(LocalizationsFinaltt));

%To save time, only investigate distances within the timescale of a
%blinking molecule

D_Overall=D_Overall(Z2<A);
Z2=Z2(Z2<A);
D_Overall=D_Overall(Z2>0);
Z2=Z2(Z2>0);

%make sure the distances are all unique, this is a silly problem, but vital

if length(unique(D_Overall))==length(D_Overall)
    good=1;
    LocalizationsFinal=LocalizationsFinaltt;
end

end

%If the true structure is known we calculate the relative score throughout.
scorey=0;
max_scorey=0;


steps=0;


if  length(TrueLocalizations)>1
    loc= True;
    fram= True_frame';
    Calculate_Score_8
    bestlik=lik+lik2;
else
    bestlik=0;
end

loc= LocalizationsFinal;
fram= Frame_Information;
Calculate_Score_8
WorstLike=lik+lik2;

Score_Newf=Inf;
%%
%Now lets determine the survival function for the distance between
%blinks

%This will be the distance cutoff for each frame separation, must be a
%monotonic decreasing function.
Surv_Blink_Dist=.5;
Surv_Blink_Distf=Surv_Blink_Dist;

Determine_Locs_w_Multi
%The other at which to link the localizations together is going to be
%important. We will permutate them with the various iterations.
Order=1:length(Frame_Information);

Orderf=Order;
Constant=0;
Constant2=0;

Traj_num=randperm(length(LocalizationsFinal),length(LocalizationsFinal));

Constantf=0;
Constant2f=0;

max_lik=-Inf;

LikHood=[];
Score=[];

Final_Localizations_Blinking_Corrected=LocalizationsFinal;
Final_Frame_Blinking_Corrected=Frame_Information;
Store_Checked_Thresholds=[];
%Cut_off_Prob2 is a threshold, where if it is larger than the prob of blink
%the localizaiton is thought to be a true localization.

scorey2=0;
Density_Calc
Timecc=0;
Timeccf=Timecc;

while steps<step_total
    
    
    Surv_Blink_Dist=Surv_Blink_Distf;
    Constant=Constantf;
    Order=Orderf;
    Timecc=Timeccf;
    
    steps=steps+1;
    
    
    vvse=randi(3);
    syc=1;
    
    if isempty(Frames_w_Multi)
        if vvse==2
            vvse=1;
        end
    end
    
    
    
    
    if vvse==1
        
        %         vvsep=randi(length(Surv_Blink_Dist));
        %         if rand<.5
        %             Surv_Blink_Dist(vvsep)=Surv_Blink_Dist(vvsep)+10*rand;
        %         else
        %             Surv_Blink_Dist(vvsep)=Surv_Blink_Dist(vvsep)-10*rand;
        %         end
        
        Surv_Blink_Dist=Surv_Blink_Distf+randn*.05;
        if Surv_Blink_Dist>1
            Surv_Blink_Dist=1;
        end
        %vvsep=randi(length(Surv_Blink_Dist));
        
        
        % Surv_Blink_Dist=Surv_Blink_Distf;
        %We start by changing 20 percent of the probabilities, much less
        %seems to just be a waste of time.
        % fracs=((steps)/(step_total+1));
        %fracs=.5;
        %
        %         if fracs<1
        %             oineppps=randperm(length(Surv_Blink_Dist), floor(length(Surv_Blink_Dist)*fracs));
        %             if isempty(oineppps)
        %                 oineppps=randi(length(Surv_Blink_Dist));
        %             end
        %         end
        %
        % Surv_Blink_Dist(oineppps)=Surv_Blink_Distf(oineppps);
        
        
        
        
        
        
    end
    
    if vvse==2
        inds42=find(Frame_Information==Frames_w_Multi(randi(length(Frames_w_Multi))));
        inds43=inds42(randperm(length(inds42)));
        Order(inds42)=Order(inds43);
    end
    
    if vvse==3
        Constant=abs(Constantf+rand*.05-.025);
    end
    
    
    %Adjust so that the order of the localizations is correct.
    %I think you need to change this to be the Localizations instead?
    Frame_Informationt=Frame_Information(Order);
    
    %Slowly adjust the fraction that is allowed to change to discourage
    %overfitting.
    
    
    [loc, fram]=Eliminate_Blinking_De_Loc14(LocalizationsFinal,  Frame_Informationt, Surv_Blink_Dist, density, Constant, Resolution, A, Deviation_in_Probability, Traj_num);
    %[loc, fram]=Threshold_multi(LocalizationsFinal, Frame_Information, Thresh_Time, Thresh_dist, density, Constants, Resolution, Constants2);
    %     if length(Surv_Blink_Dist)>length(Surv_Blink_Distf)
    %         Surv_Blink_Distf=Surv_Blink_Dist;
    %     end
    
    %Determine the Log of the liklihood for this conformation of molecules
    Calculate_Score_8
    lik+lik2-Final_lik;
    %      if  lik+lik2>=Final_lik
    %         Cut_off_Prob2=cut_t;
    %      end
    %if new best lik store everything.
    if  lik+lik2>=Final_lik || log(rand)<=(lik+lik2-Final_lik)
        
        Surv_Blink_Distf=Surv_Blink_Dist;
        Constantf=abs(Constant);
        Orderf=Order;
        Timeccf=Timecc;
        
        Final_lik=lik+lik2;
        
        
        
        if  lik+lik2>max_lik
            
            Final_Localizations_Blinking_Corrected=loc;
            Final_Frame_Blinking_Corrected=fram;
            
           
            scatter(loc(:,1),loc(:,2),10,(fram),'filled')
colormap jet
axis equal
drawnow 
            max_lik=lik+lik2;
            Final_lik=lik+lik2;
            % steps=10;
            
            if length(TrueLocalizations)>1
                [denstiy_map]=PALMplot(loc, loc1);
                scorey=sum(sum(abs(mat2gray(denstiy_map)-mat2gray(denstiy_map1)).^2));
                if scorey>max_scorey
                    [denstiy_map]=PALMplot(LocalizationsFinal, loc1);
                    scorey2=sum(sum(abs(mat2gray(denstiy_map)-mat2gray(denstiy_map1)).^2));
                    max_scorey=scorey2;
                end
                scorey2=scorey/max_scorey;
            end
            
            
            
            disp(['Fraction Complete', ' = ', num2str(Surv_Blink_Distf),' Score (If know truth)= ', num2str(scorey2),' Image= ', num2str(ksu),' MLE= ', num2str((bestlik-Final_lik))])
            %     else
            
            step_total=steps+100;
        end
        
        
        
        
        %If the true distribution of molecules is known calculate the score
        %
        
        
        %         if length(TrueLocalizations)>1
        %             [denstiy_map]=PALMplot(loc, loc1);
        %             scorey=sum(sum(abs(mat2gray(denstiy_map)-mat2gray(denstiy_map1)).^2));
        %             if scorey>max_scorey
        %                 [denstiy_map]=PALMplot(LocalizationsFinal, loc1);
        %                 scorey2=sum(sum(abs(mat2gray(denstiy_map)-mat2gray(denstiy_map1)).^2));
        %                 max_scorey=scorey2;
        %             end
        %             scorey=scorey/max_scorey;
        %         end
        %         LikHood(end+1)=lik+lik2;
        %         Score(end+1)=scorey;
        %         %disp(['Fraction Complete', ' = ', num2str(steps/step_total),' Score (If know truth)= ', num2str(scorey),' Image= ', num2str(ksu),' MLE= ', num2str((bestlik-Final_lik))])
    end
    
    if length(TrueLocalizations)>1
        [denstiy_map]=PALMplot(loc, loc1);
        scorey=sum(sum(abs(mat2gray(denstiy_map)-mat2gray(denstiy_map1)).^2));
        if scorey>max_scorey
            [denstiy_map]=PALMplot(LocalizationsFinal, loc1);
            scorey=sum(sum(abs(mat2gray(denstiy_map)-mat2gray(denstiy_map1)).^2));
            max_scorey=scorey;
        end
        scorey=scorey/max_scorey;
    end
    
    LikHood(end+1)=Final_lik;
    Score(end+1)=scorey;
    
end


end



