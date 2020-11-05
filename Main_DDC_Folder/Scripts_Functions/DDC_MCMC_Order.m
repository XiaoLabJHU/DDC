
%This is the main section of the algorithm that will eliminate blinking in
%PALM or STORM data. Please see the user guide to understand the logic of
%each section.

%Christopher Herrick Bohrer 2018, Lab of Jie Xiao and Elijah Roberts

function [Final_Localizations_Blinking_Corrected, Final_Frame_Blinking_Corrected, LikHood, Score,  Numb_of_Loc, bestlik, steps, Constantf, Orderf, Prob_dists]=DDC_MCMC(ksu, LocalizationsFinal, Frame_Information, A, Resolution, TrueLocalizations, bins, Distribution_for_Blink, steps, timer, Constantf, Orderf,  LikHood, Score, Numb_of_Loc, Prob_dists)


step_total=200;
%This is in place so that we can start off where we left off if the
%computation takes a long time.

if steps<step_total
    
    %First we organize the frames of the localizations so that they are in
    %order by frame.
    [B, Inds]=sort(Frame_Information);
    Frame_Information=Frame_Information(Inds);
    LocalizationsFinal=LocalizationsFinal(Inds,:);
    
    %We then use the following to determine the true pairwise distribution.
    DistanceMatrix={};
    for ijk=1:max(Frame_Information)+1
        DistanceMatrix{ijk}=[];
    end
    
    [DistanceMatrix]=Make_Distance_Matrix_Individual(LocalizationsFinal, DistanceMatrix, Frame_Information, max(Frame_Information));
    
    DistanceControl2=[];
    %Here we generate the True distribution using the distances that are longer
    %than the lifetime of the FP.
    for i=1:max(Frame_Information)
        if length(DistanceMatrix{i}(:))>100 && i>A
            DistanceControl2=[DistanceControl2; DistanceMatrix{i}(:)];
        end
    end
    
    
    if length(DistanceControl2)<2000
        disp('WARNING: You do not have enough data to define a true distribution, if you have split an image consider combining multiple images!!')
    end
    
    
    if length(max(DistanceControl2))==0
        error('You do not have enough localizations to define the true distribution')
    end
    
    %here I adjust the last bin so that we do not experiance a zero probability
    %in a bin, just due to noise.
    if sum(bins(2:end)>(max(DistanceControl2)))>1
        Distribution_for_Blink(bins(2:end)>(max(DistanceControl2)))=[];
        Distribution_for_Blink(end+1)=0;
        bins(bins>(max(DistanceControl2)))=[];
        bins(end+1)=Inf;
    end
    
    True_Distribuiton=((histcounts(DistanceControl2,bins,'Normalization','prob')));
    
    %Lets just convert it into probabilities
    True_Distribuiton=True_Distribuiton/sum(True_Distribuiton);
    
    %Here we will determine a rough approximation that a localization is a
    %blink, we will also generate the probability distributions that a distance
    %belonging to a blink will have, given a frame separation an certain
    %distance.
    
    %Deviation_in_Probability: is the probability that a localization is a
    %blink given a certain frame and distance away.
    
    %Prob_Distributions: is the probability of observing a distance between two localizations given that at least one of them is a blink, given the frame difference.
    
    %Dscale_store: is the proportion of the blinking distribution and of the
    %true distribution at each frame difference. With the number of blinks
    %accounted for.
    if isempty(Orderf)
        [Deviation_in_Probability, Prob_Distributions, Dscale_store]=Determine_Deviation_in_Probability8(True_Distribuiton,Distribution_for_Blink, bins, A, LocalizationsFinal, Frame_Information, [], [],Resolution);
        Prob_Distributions=log(Prob_Distributions);
        Prob_dists.Deviation_in_Probability=Deviation_in_Probability;
        Prob_dists.Prob_Distributions=Prob_Distributions;
        Prob_dists.Dscale_store=Dscale_store;
    else
        Deviation_in_Probability=Prob_dists.Deviation_in_Probability;
        Prob_Distributions=Prob_dists.Prob_Distributions;
        Dscale_store=Prob_dists.Dscale_store;
    end
    
    %True localizations if availible from simulaiton.
    if length(TrueLocalizations)>1
        True=[];
        True_frame=[];
        
        vvseo=zeros(1,length(LocalizationsFinal(:,1)));
        for kksd=1:length(LocalizationsFinal(:,1))
            
            if sum(find(LocalizationsFinal(kksd,1)==TrueLocalizations(:,1)))>0
                True=[True; LocalizationsFinal(kksd,:)];
                True_frame=[True_frame; Frame_Information(kksd)];
                vvseo(kksd)=1;
            else
                
            end
        end
        
        loc1=True;
        loc=True;
        %Generate the true image to calculate the relative score if simulation.
        [denstiy_map1]=PALMplot(loc, loc1);
    end
    
    
    %here we start calculating the pairwise distances.
    D_Overall = (pdist(LocalizationsFinal));
    
    %Calculate the pairwise distances and the overall frame differences
    Z2 = pdist([[1:length(Frame_Information)]'*0,Frame_Information(:)]);
    
    D_Overall=D_Overall(Z2<A);
    Z2=Z2(Z2<A);
    D_Overall=D_Overall(Z2>0);
    Z2=Z2(Z2>0);
    
    %Calculate the likilihood of the true localizations if it is availible.
    if  length(TrueLocalizations)>1
        loc= True;
        fram= True_frame';
        Calculate_Score_8
        bestlik=lik+lik2;
    else
        bestlik=0;
    end
    
    Z2 = pdist([[1:length(Frame_Information)]'*0,Frame_Information(:)]);
    
    % for a stupid reason some people round their localizations and this can screw things up, here I am correcting for this issue.
    LocalizationsFinaltt=LocalizationsFinal;
    good=0;
    while good==0
        
        LocalizationsFinaltt=LocalizationsFinal+rand(size(LocalizationsFinal))*.01;
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
    
    if length(TrueLocalizations)>1
        [denstiy_map]=PALMplot(LocalizationsFinal, loc1);
        scorey=sum(sum(abs(mat2gray(denstiy_map)-mat2gray(denstiy_map1)).^2));
        max_scorey=scorey;
    end
    
    
    
    Final_lik=-Inf;
    lik=999999;
    
    loc= LocalizationsFinal;
    fram= Frame_Information;
    Calculate_Score_8
    WorstLike=lik+lik2;
    
    %%
    %Here we determine if there is any problem of two localizations having the
    %same frame and within resolution of each other. see the following script
    %for more details
    
    
    %The other at which to link the localizations together is going to be
    %important. We will permutate them with the various iterations.
    % Determine_Locs_w_Multi2
    
    %Here we determine the relative density of the localizations using the
    %localizations that are separated by at least the lifetime of the FP.
    %Density_Calc
    
    unsorted = 1:length(LocalizationsFinal);
    
    Traj_num=randperm(length(LocalizationsFinal),length(LocalizationsFinal));
    
    
    if isempty(Orderf)
        Order=1:length(Frame_Information);
        
        Orderf=Order;
        %Here you define the number of terms that you are going to use for the
        %arbitrary polynomial.
        Constant=(1:length(Frame_Information))*0;
        Constantf=Constant;
        
        Numb_of_Loc=[];
        LikHood=[];
        Score=[];
        max_lik=-Inf;
        
        Final_Localizations_Blinking_Corrected=LocalizationsFinal;
        Final_Frame_Blinking_Corrected=Frame_Information;
    else
        max_lik=max(LikHood);
        
        %Adjust so that the order of the localizations is correct.
        LocalizationsFinalt=LocalizationsFinal(Orderf,:);
        [Final_Localizations_Blinking_Corrected, Final_Frame_Blinking_Corrected]=Eliminate_Blinking_De_Loc15(LocalizationsFinalt,  Frame_Information, Constantf(Orderf), Resolution, A, Deviation_in_Probability, Traj_num);
        
    end
    
    
    
    tic
    
    
    while steps<step_total
        
        %if the time for this one is going to long stop and output the data
        %This is in place if the computation could take a long time.
        if toc>timer
            disp('Stopping to save some of the data')
            break
        end
        
        Constant=Constantf;
        Order=Orderf;
        steps=steps+1;
        
        %Here we go through and adjust constant and the order in which we
        %link the localizations into trajectories and make it so that we
        %have to investigate a new conforamtion of localizaitons.
        
        loc2=loc;
        countss=0;
        while isequal(loc, loc2)
            
            inds42=randi([2,length(Frame_Information)-1],1);
            if rand<.5
                inds43=inds42+1;
            else
                inds43=inds42-1;
            end
            Ordert=Order;
            Order(inds42)=Ordert(inds43);
            Order(inds43)=Ordert(inds42);
            
            
            %Adjust so that the order of the localizations is correct.
            LocalizationsFinalt=LocalizationsFinal(Order,:);
            [loc, fram]=Eliminate_Blinking_De_Loc15(LocalizationsFinalt,  Frame_Information(Order), Constant(Order)*0, Resolution, A, Deviation_in_Probability, Traj_num);
            
        end
        
        %Determine the Log of the liklihood for this conformation of molecules
        Calculate_Score_8
        
        %Here is the heart of the MCMC approach.
        if  lik+lik2>=Final_lik || log(rand)<=(lik+lik2-max_lik)
            
            Constantf=(Constant);
            Orderf=Order;
            Final_lik=lik+lik2;
            
            %Store the conformation of localizations if it is the max lik.
            if  lik+lik2>max_lik
                
                Final_Localizations_Blinking_Corrected=loc;
                Final_Frame_Blinking_Corrected=fram;
                
                %                 %If you are not in // and want to look at the
                %                 localizations in real time.
                %                 scatter(loc(:,1),loc(:,2),10,(fram),'filled')
                %                 colormap jet
                %                 axis equal
                %                 drawnow
                
                max_lik=lik+lik2;
                Final_lik=lik+lik2;
                steps=10;
            end
            
            
            %Here we calculate the relative score if you are analyzing
            %simulation data.
            if length(TrueLocalizations)>1
                [denstiy_map]=PALMplot(loc, loc1);
                scorey11=sum(sum(abs(mat2gray(denstiy_map)-mat2gray(denstiy_map1)).^2));
                scorey=scorey11/max_scorey;
            end
            
            LikHood(end+1)=Final_lik;
            Score(end+1)=scorey;
            Numb_of_Loc(end+1)=length(loc);
        end
        
        disp(['Fraction Complete', ' = ', num2str(steps/step_total),' Score (If know truth)= ', num2str(scorey),' Image= ', num2str(ksu),' MLE= ', num2str((bestlik-Final_lik))])
        
    end
    
    
end


end
