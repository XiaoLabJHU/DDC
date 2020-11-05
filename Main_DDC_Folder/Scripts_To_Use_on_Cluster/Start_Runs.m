%% DDC (i.e. The Blinking Algorithm)
% Christopher H. Bohrer 2019
% If you have any questions please contact me at cbohrer1@jhmi.edu
% Really just send me an email!!

% This is a program which will correct for the blinking of a FP or dye
% without the use of any inputs other than the data iself. Well it does
% need two, the frame difference at which the true pairwise distributions
% are obtained and the resolution of the SMLM.

% If you used different imaging
% conditions or did anything different accross the different images they
% should be imaged separatly.

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!Make sure the localizations are in nm!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Please make sure the data is stored in a cell with the name
% LocalizationsFinal and the frame information in the cell Frame_Information

% See the USER GUIDE for an example on how to use this code.
clear
files=dir('Split*');
% First clear the workspace.
for oinepppq=1:length(files)
    TrueLocalizations=[];
    files=dir('Split*');
    filename=files(oinepppq).name
    load(files(oinepppq).name)
    
    tryers=dir(['A*',filename,'*']);
    if length(tryers)==0
        try
            % Set to 1 if using a cluster to analyze the data. You will need to modify
            % the script 'Run_cluster.m' depending upon your particular case.
            cluster=0;
            
            % If you want to investigate how simulated data performs you can look at the
            % relative score. Otherwise we create a blank array, as the only way you can
            % know the true localizations is with simulation data.
            
           
            addonarray=[];
            addonarray(500)=0;
            
            % This will adjust your path based off of whether you are using a cluster or
            % not. You may need to adjust the names of some strings in the following if
            % you are using a cluster. It is very intuative, if not just leave it alone
            %Run_cluster
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Resolution=20;% Determined using Determine_Res.m
            stepper=150; % Maximum number of MCMC steps
            N_f=700; % Determined using Determine_Res.m
            
            
            
            % This is only here if you DO NOT create simulation data.
            if isempty(TrueLocalizations)
                TrueLocalizations={};
                for ijk=1:length(Frame_Information)
                    TrueLocalizations{ijk}=[];
                    
                end
            end
            
            % If the localizations are in 2D adjust the analysis
            for sdfv=1:length(LocalizationsFinal)
                if min(size(LocalizationsFinal{sdfv}))<3
                    LocalizationsFinal{sdfv}(:,3)=LocalizationsFinal{sdfv}(:,2)*0;
                end
            end
            
            % The program will save out the data after analysis.
            % It will have the following name and date added to the name.
            String=['Analyzed_Time',datestr((datetime('today'))),'_',filename];
            
            %These arrays will save the results! They will provide the final cordinates
            %of the images.
            Final_Loc_Blinking_Corrected={};%This is the set that contains the localizations that are considered the real localizations
            Final_Frame_Blinking_Corrected={};%and their frames
            
            %% Determining the Blinking Distribution
            % This is the first step of DDC, where the blinking distribution is actually
            % determined. This step will determine P_{blink} (See Supporting material)
            
            % Here we determine the probability distribution that a molecule will blink
            % a certain distance away. It is calculated by compairing the true
            % distributions from all cells, with the true distributions from all cells.
            
            % Note here all of the images within the LocalizationFinal structure are
            % supposed to have the same resolution. If you used different imaging
            % conditions or did anything different accross the different images they
            % should be imaged separatly.
            
            [bins, Distribution_for_Blink, ~, Resolution, X_overall, M_mat] =...
                Determine_Blinking_Distribution5(LocalizationsFinal,...
                Frame_Information, N_f, Resolution);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            timer=3600*.1;%The amount of time to spend investigating each image before saving output.
            
            % Here is some information that is used by DDC. You do not need to worry
            % about the specifics of these cells. Though, if you do want to understand
            % them you can see the USER GUIDE!
            
            Constantf={}; % This is known as kappa(density) within the supporting material
            Constantf2={}; % This is known as kappa2(frame) within the supporting material
            Orderf={}; % This is the order of the localizations in which they are linked into trajectories
            Step=[]; % Stores where each image is within the MCMC phase space search
            Lik={}; % The likelihood at that step
            RelScore={}; % The image error if the true localizations are known (assuming 40nm resolution for the SMLM data)
            Numb_of_Loc={}; % The number of localizations at each step of the MCMC approach for each image
            Prob_dists={}; % The probability distributions for each image used to calculate the likelihood, so that we do not need to calculate them each time.
            
            for ijk=1:length(Frame_Information)
                Orderf{ijk}=[];
                Lik{ijk}=[];
                RelScore{ijk}=[];
                Numb_of_Loc{ijk}=[];
                Constantf{ijk}=[];
                Constantf2{ijk}=[];
                Step(ijk)=0;
                
                Prob_dists{ijk}.Deviation_in_Probability=[];
                Prob_dists{ijk}.Prob_Distributions=[];
                Prob_dists{ijk}.Dscale_store=[];
            end
            
            save([String])
        end
    else
        disp('Allready did that')
    end
end
%% Blinking Elimination
% This is the main section of DDC, where the blinking is actually
% eliminated.



%{
while sum(Step<stepper)>0
    % In this part of the code we will eliminate blinking in individual images.
    hey=1
    parfor ksu=1:length(LocalizationsFinal)
        
        % This is to guarantee that there are not to many localizations used
        % during the analysis with DDC, if you have more than 6000
        % localizations per image DDC's effectivness will be greatly
        % depleted;
        
        
        % If it is greater than 7000 do not analyze it. Consider breaking up
        % the images.
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
            
            % Store out the imformation from the MCMC approach.
            Constantf{ksu}=Constant;
            Constantf2{ksu}=Constant2;
            Orderf{ksu}=Order;
            Step(ksu)=steps;
            Lik{ksu}=LikHood;
            RelScore{ksu}=Score;
            Numb_of_Loc{ksu}=Numb;
            
            % If the best lik was reached store the new set
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

clear

%}