%% Determine_Res.m
%
%This is a script that allows you to determine the Bin Resolution for
% application of DDC (i.e. when the bin width is sufficient to result in
% a good approximation of the true pairwise distance distribution).
%
% See 'Determining the Bin Resolution" section of the USER GUIDE.

%usually the resol
clear
Res=80;  % Define the max number of frames to probe
N_f=200; % This is N as determined in "Determine_N.m"

% We will only investigate the first 10% of the images. if you want more,
% you can change it. See "Determine the max size (D_maxi) of the 
% image being processed" section of the code below.
image_frac = 0.1;

% Hit run and we select the data to load (i.e. Split_*.m)
% Add the proper folders to the path


Location=pwd;
addpath([Location, '/Scripts_Functions'])

% This will allow you to pick out your own file
[filename, pathname] = uigetfile({'*.mat'}, 'Select HMM .mat file');
if ( filename == 0 )
    disp('Error! No (or wrong) file selected!')
    return
end

% This will load in all of the localizations from your structure file
full_filename = [ pathname, filename ];
load(full_filename,'Frame_Information','filename','TrueLocalizations','LocalizationsFinal');

%% Determine the max size (D_maxi) of the image being processed

D_maxi=0;
for i=1:ceil(length(LocalizationsFinal))
    
    D = (pdist(LocalizationsFinal{i}));
    if D_maxi<max(D)
        D_maxi=max(D);
    end
    
end

bins=[0:Res:D_maxi-Res, Inf];
Res_Store=[];
for i=1:ceil(length(LocalizationsFinal)*image_frac)
    i/ceil(length(LocalizationsFinal))
    Z2 = pdist([[1:length(Frame_Information{i})]'*0,Frame_Information{i}(:)]);
    D = (pdist(LocalizationsFinal{i}));
    
    D=D((Z2)>N_f);
    Res_Store=[Res_Store; histcounts(D,bins, 'Normalization','Prob')];
    
end

%% Ploting TPDD using All Localizations

% See Figure 3 (top) of the USER GUIDE.

subplot(2,1,1)
plot((Res_Store'))
xlabel('Distance','FontSize',20)
ylabel('Prob','FontSize',20)
title('Using All Localizations','FontSize',30)
axis tight
legend

%% Groud Truth Comparison

% This section compares simulated distributions to their ground truth
% counterparts. Ground truth localizations should be stored in a
% TrueLocalizations cell (see "Determining the Bin Resolution" section of
% the USER GUIDE for more information.

if length(TrueLocalizations)>0
    
    Res_Store2=[];
    
    
    for i=1:ceil(length(LocalizationsFinal)*image_frac)
        True=[];
        True_frame=[];
        
        vvseo=zeros(1,length(LocalizationsFinal{i}(:,1)));
        for kksd=1:length(LocalizationsFinal{i}(:,1))
            
            if sum(find(LocalizationsFinal{i}(kksd,1)==TrueLocalizations{i}(:,1)))>0
                True=[True; LocalizationsFinal{i}(kksd,:)];
                True_frame=[True_frame; Frame_Information{i}(kksd)];
            else
                
            end
        end
        
        i/ceil(length(LocalizationsFinal))
        Z2 = pdist([[1:length(True_frame)]'*0,True_frame]);
        D = (pdist(True));
        
        D=D((Z2)>N_f);
        Res_Store2=[Res_Store2; histcounts(D,bins, 'Normalization','Prob')];
        
    end
    
%% Plotting TPDD using True Localiations
    
    % See Figure 3 (bottom) in the USER GUIDE.
    
    subplot(2,1,2)
    plot((Res_Store2'))
    axis tight
    xlabel('Distance','FontSize',20)
    ylabel('Prob','FontSize',20)
    title('Using True Localizations','FontSize',30)
    legend
end
