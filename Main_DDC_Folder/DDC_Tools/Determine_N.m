%% Determine_N.m
%
% This will be a script that will allow you to determine "N" (i.e. when the
% pairwise distance distributions reach a steady state) resulting in the
% true pairwise distribution.
%
% See "Determining the Frame Difference (N)" section of the USER GUIDE.


clear

% Define NFTP (maximum number of iternations to probe)

NFTP=1000;
Gap=10;%This is how far each frame difference jump will be, the larger this number 
%the faster this script will go, but the less clear the final graph will be


% Hit run and select the data to load (i.e. Split_*.m)
% Add the proper folders to the path


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
%load(full_filename,'Frame_Information','TrueLocalizations','LocalizationsFinal');
load(full_filename,'Frame_Information','filename','TrueLocalizations','LocalizationsFinal');

%% Calculation of Cumulative Distribution Function

Cum_Sum_Store=[];
means=[];
Framestore=[];
bins=[0:250:5000, Inf];
disp('Determining N')
for iis=1:Gap:NFTP
    disp(['Fraction Done=',num2str(iis/NFTP)])
    Total_Blink=[];
    for i=1:ceil(length(LocalizationsFinal))
        if length(LocalizationsFinal{i})<6000
        % Gather the distance distributions from 10% of the images to speed up
        % the analysis. This should be sufficient to to define Z.
        
        %if length(LocalizationsFinal{i})>2000
        Z2 = pdist([[1:length(Frame_Information{i})]'*0,Frame_Information{i}(:)]);
        D = (pdist(LocalizationsFinal{i}));
        D = D(abs(Z2-iis)==0);
        
        Total_Blink=[Total_Blink, D];
        % end
        end
    end
    Cum_Sum_Store=[Cum_Sum_Store; histcounts(Total_Blink,bins,'Normalization','cdf')];
    Framestore(end+1)=iis;
end

%% Calculation of Z (as in main text)

Z=[];

for iis2=1:size(Cum_Sum_Store,1)
    iis2
    Z(end+1)=sum(abs(Cum_Sum_Store(1,:)-Cum_Sum_Store(iis2,:)));
end

%% Plotting Z vs. Frame Difference

% See Figure 2 in the USER GUIDE

figure(2)
plot(Framestore,Z,'LineWidth',1.25)
xlabel('Frame')
ylabel('Z')
set(gca,'FontSize',20)
save temp_file_with_Z_vs_frame
