
%This is the main section of the algorithm that will eliminate blinking in
%PALM or STORM data. Please see the user guide to understand the logic of
%each section/ Supporting material.

%Christopher Herrick Bohrer 2018, Lab of Jie Xiao and Elijah Roberts

function [Traj2, Final_Localizations_Blinking_Corrected, Final_Frame_Blinking_Corrected, LikHood, Score,  Numb_of_Loc, bestlik, steps, Constantf, Orderf, Prob_dists, Constant2ft]=DDC_MCMC(ksu, LocalizationsFinal, Frame_Information, A, Resolution, TrueLocalizations, bins, Distribution_for_Blink, steps, timer, Constantf, Orderf,  LikHood, Score, Numb_of_Loc, Prob_dists, step_total, addon, Constant2ft, X_overall, M_mat)

%Load in the probability distribuitons used in the likelihood calculations
%and the matrix M, to approx the probability that a loc is a blink.
Prob_dists.Deviation_in_Probability=full(Prob_dists.Deviation_in_Probability);
Prob_dists.Prob_Distributions=full(Prob_dists.Prob_Distributions);


if steps<step_total
    
    
    %First we organize the frames of the localizations so that they are in
    %order by frame.
    [~, Inds]=sort(Frame_Information);
    Frame_Information=Frame_Information(Inds);
    LocalizationsFinal=LocalizationsFinal(Inds,:);
    
    %We then use the following to determine the true pairwise distribution.
    DistanceMatrix={};
    for ijk=1:max(Frame_Information)+1
        DistanceMatrix{ijk}=[];
    end
    
    %This will go through and grab the distances for the real pairwise
    %distance distribution.
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
    
    
%     
%     if isempty(max(DistanceControl2))
%         disp('You do not have enough localizations to define the true distribution')
%         for i=1:max(Frame_Information)
%             if length(DistanceMatrix{i}(:))>1 && i>A
%                 DistanceControl2=[DistanceControl2; DistanceMatrix{i}(:)];
%             end
%         end
%     end
%     
    
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
    
    
    %By using larger bins for the fitting of X (See Supporting Material) you
    %obtain less error.
    
    Distribution_for_Blink22=[];
    for oine=1:2:length(Distribution_for_Blink)
        if length(Distribution_for_Blink)>=oine+1
            Distribution_for_Blink22(end+1)=(Distribution_for_Blink(oine)+Distribution_for_Blink(oine+1));%+Distribution_for_Blink(oine+2)+Distribution_for_Blink(oine+3));
        end
    end
    
    
    bins22=[];
    for oine=1:2:length(bins)
        if length(Distribution_for_Blink)>=oine+1
            bins22(end+1)=bins(oine);
        end
    end
    
    
    bins22(end+1)=Inf;
    
    True_Distribuiton22=((histcounts(DistanceControl2,bins22,'Normalization','prob')));
    
    %Lets just convert it into probabilities
    True_Distribuiton22=True_Distribuiton22/sum(True_Distribuiton22);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %Here we will determine a rough approximation that a localization is a
    %blink, we will also generate the probability distributions that a distance
    %belonging to a blink will have, given a frame separation an certain
    %distance.
    
    %Deviation_in_Probability: is the probability that a localization is a
    %blink given a certain frame and distance away. The matrix M
    
    %Prob_Distributions: is the probability of observing a distance between
    %two localizations given that at least one of them is a blink, given
    %the frame difference. This is the probability distribution P_B1 within
    %the main text.
    
    %Dscale_store: is the proportion of the blinking distribution and of the
    %true distribution at each frame difference. With the number of blinks
    %accounted for. This is X, from the fitting in the supporting material of 
    %the main text.
    
    
    if isempty(Orderf)
        [Deviation_in_Probability, Prob_Distributions, Dscale_store]=Determine_Deviation_in_Probability8(True_Distribuiton, Distribution_for_Blink, bins, A, LocalizationsFinal, Frame_Information, [], [],Resolution,True_Distribuiton22,Distribution_for_Blink22, bins22, X_overall, M_mat);
        Prob_Distributions=log(Prob_Distributions);
        Prob_dists.Deviation_in_Probability=Deviation_in_Probability;
        Prob_dists.Prob_Distributions=Prob_Distributions;
        Prob_dists.Dscale_store=Dscale_store;
    else
        Deviation_in_Probability=Prob_dists.Deviation_in_Probability;
        Prob_Distributions=Prob_dists.Prob_Distributions;
        Dscale_store=Prob_dists.Dscale_store;
    end
    
    
    
    %True localizations if availible from simulaiton, grabe frames of the
    %true localizations.
    if length(TrueLocalizations)>1
        
        True=[];
        True_frame=[];
        vvseo=zeros(1,length(LocalizationsFinal(:,1)));
        
        for kksd=1:length(LocalizationsFinal(:,1))
            
            if sum(find(LocalizationsFinal(kksd,1)==TrueLocalizations(:,1)))>0
                
                True=[True; LocalizationsFinal(kksd,:)];
                True_frame=[True_frame; Frame_Information(kksd)];
                vvseo(kksd)=1;
                
            end
        end
        
        loc1=True;
        loc=True;
        
        %Generate the true image to calculate the relative score if
        %simulation. This is assuming a resolution of 40 nm. You may want
        %to be more specific.
        [denstiy_map1]=PALMplot(loc(:,1:2), loc1(:,1:2));
        
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
    
    % for some reason some people round their localizations and this can screw things up, here I am correcting for this issue.

    
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
        [denstiy_map]=PALMplot(LocalizationsFinal(:,1:2), loc1(:,1:2));
        %denstiy_map=
        % denstiy_map=denstiy_map(mmask>0);
        scorey=sum(sum(abs(mat2gray(denstiy_map)-mat2gray(denstiy_map1)).^2));
        max_scorey=scorey;
    end
    
    
    
    Final_lik=-Inf;
    lik=999999;
    
    
    loc= LocalizationsFinal;
    fram= Frame_Information;
    %Calculate_Score_8
   % WorstLike=lik+lik2;
    
    %%
    %Here we determine if there is any problem of two localizations having the
    %same frame and within resolution of each other. see the following script
    %for more details
    
    
    %The other at which to link the localizations together is going to be
    %important. We will permutate them with the various iterations.
    Determine_Locs_w_Multi2
    
    %Here we determine the relative density of the localizations using the
    %localizations that are separated by at least the lifetime of the FP.
    Density_Calc
    
    unsorted = 1:length(LocalizationsFinal);
    
    %To store the trajectory of the individual localization (make sure you
    %check this is 
    Traj_num=randperm(length(LocalizationsFinal),length(LocalizationsFinal));
    Traj=Traj_num;
    Traj2=Traj;
    
    if isempty(Orderf)
        
        %If this is the first time running DDC on this data, we are going
        %to have to make some arrays and stuff for the MCMC approach. 
        
        Order=1:length(Frame_Information); 
        Orderf=Order;%Whenever there is an f, that means that it is the final or the stored one that resulted from maximizing the lik
        Constant=(1:length(Frame_Information))*0;%Initial starting conditions for kappas are equal to zero
        Constantf=Constant;
        Constant2=Constant;
        Constant2ft=Constant;
        Numb_of_Loc=[];
        LikHood=[];
        Score=[];
        max_lik=-Inf;
        Final_Localizations_Blinking_Corrected=LocalizationsFinal;
        Final_Frame_Blinking_Corrected=Frame_Information;
    else
        
        %If the data is just continuing to run, then we have to figure out
        %where the MCMC approach left off.
        max_lik=max(LikHood);
        %Adjust so that the order of the localizations is correct.
        LocalizationsFinalt=LocalizationsFinal(Orderf,:);
        kappas=Constantf(Orderf)+Constant2ft(Orderf);
        kappas(kappas<-1)=-.9999;
        [Final_Localizations_Blinking_Corrected, Final_Frame_Blinking_Corrected, Traj2]=Eliminate_Blinking_De_Loc15(LocalizationsFinalt,  Frame_Information, kappas, Resolution, A, Deviation_in_Probability, Traj_num);
        Traj2=Traj2(Orderf);
    end
    
    
    
    tic
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Here is the actual MCMC phase space search. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    while steps<step_total
        
        %if the time for this one is going to long stop and output the data
        %This is in place if the computation could take a long time.
        if toc>timer
            disp('Stopping to save some of the data')
            break
        end
        
        Constant=Constantf;
        Constant2=Constant2ft;
        Order=Orderf;
        steps=steps+1;
        
        
        %Here we go through and adjust constant ( (kappa(density) and kappa2(frame) ) and the order in which we
        %link the localizations into trajectories and make it so that we
        %have to investigate a new conforamtion of localizaitons.
        
        loc2=loc;
      
        %make sure that we for each step we at least investigate a new
        %conformation of localizations!! 
        
        while isequal(loc, loc2)
            
            %generate a random interager, if 1 adjust kappa(density) if 2
            %adjust kappa2(frame) if 3 adjust order. 
            vvse=randi(3);
            if isempty(Frames_w_Multi)
                if vvse==2
                    vvse=1;
                end
            end
            
            %This code also make sure that the functions (kappa(density) kappa(frame) are monotonically
            %increasing.
            
            %This will adjust for the density of the image kappa(density).
            if vvse==1
                
                %Again if we are adjusting the constant monotonic function
                %we make sure that it is different at the end of this step.
                
                if steps>1
                    
                    [~,III2] = sort( density);
                    III(III2) = unsorted;
                    Constantft=Constant;
                    
                    %Make sure that we investigate a new conformation of
                    %the function 
                    while isequal(Constant,Constantf) || isequal(Constant,Constantft)
                        
                        %Adjust random position of the function 
                        va=randi(length(Constant));
                        Constant(va)=Constant(va)+randn*.1;
                        
                        %Ensure that it is a monotonic function, can do
                        %this in two different ways, therefore choose which
                        %way randomly.
                        if rand<.5
                            for iserfv =2: numel(Constant)
                                if Constant(III==(iserfv - 1)) > Constant(III==iserfv)
                                    Constant(III==(iserfv))= Constant(III==(iserfv-1));
                                end
                            end
                        else
                            for iserfv = numel(Constant):-1:2
                                if Constant(III==(iserfv - 1)) > Constant(III==iserfv)
                                    Constant(III==(iserfv-1))= Constant(III==(iserfv));
                                end
                            end
                        end
                    end
                    
                    %Set the constraint that it cannot go below -.9, as
                    %this limits the phase space search and if there is a
                    %negative number in the numerator within the script
                    %Eliminat_De_le then there will be issues.
                    Constant(Constant<-.9)=-.9;
                    
                    %Add one more constriant, that the mean over all values
                    %must equal zero. This maintains that the phase space
                    %search will not deviate to far from the more probable
                    %conformations of localizations. 
                    if mean(Constant)<0
                        Constant=Constant+abs(mean(Constant));
                    else
                        Constant=Constant-abs(mean(Constant));
                    end
                end
            end

            
            if  vvse==3
                %Again if we are adjusting the Constant2 monotonic function
                %we make sure that it is different at the end of this step.
                if steps>1
                    
                    Constant2ft=Constant2;
                    
                    while isequal(Constant2,Constant2ft) || isequal(Constant2,Constant2ft)

                        va=randi(length(Constant2));
                        Constant2(va)=Constant2(va)+randn*.1;

                        if rand<.5
                            for iserfv =2: numel(Constant2)
                                if Constant2((iserfv - 1)) < Constant2(iserfv)
                                    Constant2((iserfv))= Constant2((iserfv-1));
                                end
                            end
                        else
                            for iserfv = numel(Constant2):-1:2
                                if Constant2((iserfv - 1)) < Constant2(iserfv)
                                    Constant2((iserfv-1))= Constant2((iserfv));
                                end
                            end
                        end
                        
                        Constant2(Constant2<-.9)=-.9;
                        
                        if mean(Constant2)<0
                            Constant2=Constant2+abs(mean(Constant2));
                        else
                            Constant2=Constant2-abs(mean(Constant2));
                        end
                    end
                end
            end
            
            kappas=Constant+Constant2;
            kappas(kappas<-1)=-.9999;
            
            
            %Here we have the option where DDC can adjust the order in
            %which the localizations are linked. 
            if vvse==2
                
                inds42=find(Frame_Information==Frames_w_Multi(randi(length(Frames_w_Multi))));
                inds43=inds42(randperm(length(inds42)));
                Order(inds42)=Order(inds43);
               
            end
            
            %Adjust so that the order of the localizations is correct.
            LocalizationsFinalt=LocalizationsFinal(Order,:);
            
            
            %This is alg 1, linking localizations into trajectories from
            %the supporting material!! 
            [loc, fram, Traj]=Eliminate_Blinking_De_Loc15(LocalizationsFinalt,  Frame_Information, kappas(Order), Resolution, A, Deviation_in_Probability, Traj_num);
            Traj=Traj(Order);
            
        end
        
        %Determine the Log of the liklihood for this conformation of molecules
        Calculate_Score_8
        
        %Here is the heart of the MCMC approach. Store the parameters if
        %they meet this criteria
        if  lik+lik2>=Final_lik || log(rand)<=(lik+lik2-Final_lik)
            
            Constantf=(Constant);
            Constant2ft=(Constant2);
            Orderf=Order;
            Final_lik=lik+lik2;
            
            %Store the conformation of localizations if it is the max lik.
            if  lik+lik2>max_lik
                
                Final_Localizations_Blinking_Corrected=loc;
                Final_Frame_Blinking_Corrected=fram;
                Traj2(Inds)=Traj;
                %If you are not in // and want to look at the
                % localizations in real time.
                scatter(loc(:,1),loc(:,2),10,(fram),'filled')
                colormap jet
                axis equal
                drawnow
                
                max_lik=lik+lik2;
                Final_lik=lik+lik2;
                steps=10;
            end
            
            
            %Here we calculate the relative score if you are analyzing
            %simulation data.
            if length(TrueLocalizations)>1
                [denstiy_map]=PALMplot(loc(:,1:2), loc1(:,1:2));
                %denstiy_map=denstiy_map(mmask>0);
                scorey11=sum(sum(abs(mat2gray(denstiy_map)-mat2gray(denstiy_map1)).^2));
                scorey=scorey11/max_scorey;
            end
            
            LikHood(end+1)=Final_lik;
            Score(end+1)=scorey;
            Numb_of_Loc(end+1)=length(loc);
        end
        
       %disp(['Fraction Complete', ' = ', num2str(steps/step_total),' Score (If know truth)= ', num2str(scorey),' Image= ', num2str(ksu),' MLE= ', num2str((bestlik-Final_lik))])
        disp(['Fraction Complete', ' = ', num2str(steps/step_total),' Score (If know truth)= ', num2str(scorey),' Image= ', num2str(ksu)])
    end
    
    
end


end
