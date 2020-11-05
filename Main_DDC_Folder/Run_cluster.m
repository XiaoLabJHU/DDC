if cluster==1;
    
    Location=pwd;
    parpool(number_of_workers)
    
    %If you are using a cluster you will have to define your filename.
    
    filename=['long_time_2_16000.mat'];
    pathname=Location;
    
    %Add global scripts to path.
    addpath([Location, '/Scripts_Functions'])
    
    %This will load in all of the localizations from your structure file
    full_filename = [ pathname, '/', filename ];
    %load(full_filename, 'Frame_Information','filename','TrueLocalizations','LocalizationsFinal');
    load(full_filename);
else
    
    %Add the proper folders to the path
    Location=pwd;
    addpath([Location, '/Scripts_Functions'])
    
    %This will allow you to pick out your own file
    [filename, pathname] = uigetfile({'*.mat'}, 'Select HMM .mat file');
    if ( filename == 0 )
        disp('Error! No (or wrong) file selected!')
        return
    end
    
    %This will load in all of the localizations from your structure file
    full_filename = [ pathname, filename ];
    %load(full_filename,'Frame_Information','filename','TrueLocalizations','LocalizationsFinal');
    load(full_filename);
end