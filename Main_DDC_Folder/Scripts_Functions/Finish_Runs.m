clear
%This is going to be a script to finish runs
for oinepppq=1:1000
    clear
    files=dir('An*');
    load(files(randi(length(files))).name)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is the main section of DDC, where the blinking is actually eliminated

if sum(Step<stepper)>0
    %In this part of the code we will eliminate blinking in individual images.
    
    parfor ksu=1:length(LocalizationsFinal)
        
        %This is to guarantee that there are not to many localizations used
        %during the analysis with DDC, if you have more than 6000
        %localizations per image DDC's effectivness will be greatly
        %depleted; 
      
        
        %If it is greater than 7000 do not analyze it. Consider breaking up
        %the images. 
        if length(LocalizationsFinal{ksu})>7000
            disp('Warning: you may have to many localizations for each image, consider splitting it up')
        end
        
        if length(LocalizationsFinal{ksu})>10 && Step(ksu)<stepper
        
            [Final_Localizations_Blinking_Corrected2, Final_Frame_Blinking_Corrected2, LikHood, Score,  Numb, bestlik, steps, Constant, Order, Prob_dists2, Constant2]=DDC_MCMC(ksu, LocalizationsFinal{ksu}, round(Frame_Information{ksu}), N_f, Resolution, TrueLocalizations{ksu}, bins, Distribution_for_Blink, Step(ksu), timer, Constantf{ksu}, Orderf{ksu}, Lik{ksu}, RelScore{ksu}, Numb_of_Loc{ksu}, Prob_dists{ksu}, stepper, addonarray(ksu), Constantf2{ksu}, X_overall, M_mat);
            
            maxtemp=max(Lik{ksu});
            if isempty(maxtemp)
                maxtemp=-Inf;
            end
            
            Prob_dists{ksu}.Deviation_in_Probability=sparse(Prob_dists2.Deviation_in_Probability); %Note here that Deviation_in_Probability is actually M matrix from the supporting material
            Prob_dists{ksu}.Prob_Distributions=sparse(Prob_dists2.Prob_Distributions); %These distributions are the P_B1 from the supporting material given the frame difference of the row
            Prob_dists{ksu}.Dscale_store=Prob_dists2.Dscale_store;% this is the proportion of P_R and P_{blink}, aka X from the supporting material. 
            
            %Store out the imformation from the MCMC approach. 
            Constantf{ksu}=Constant;
            Constantf2{ksu}=Constant2;
            Orderf{ksu}=Order;
            Step(ksu)=steps;
            Lik{ksu}=LikHood;
            RelScore{ksu}=Score;
            Numb_of_Loc{ksu}=Numb;
            
            %If the best lik was reached store the new set
            if maxtemp<max(LikHood)
                Final_Loc_Blinking_Corrected{ksu}=Final_Localizations_Blinking_Corrected2;
                Final_Frame_Blinking_Corrected{ksu}=Final_Frame_Blinking_Corrected2;
            end
            
            
        else
            Step(ksu)=stepper;
        end
    
    end
    save([String])
    
end




end