%% Combine_Images_Multiple_Images.m

%This script recombines your images after DDC has been run

% The target number of localizations is  min_loc and some buffer will be
% added to each side of the images to avoid boundary effects.

clear
% Select the image (i.e. Analyzed_*.mat file) you want to combine.
% Pick your own file
[filename, pathname] = uigetfile({'*.mat'}, 'Select .mat file');
if ( filename == 0 )
    disp('Error! No (or wrong) file selected!')
    return
end

Cutsstore=[];
% Load in all of the localizations from your structure file
full_filename = [ pathname, filename ];
load(full_filename);

Condition=filename;

 %Final_Loc_Blinking_Corrected=LocalizationsFinal;
 %Final_Frame_Blinking_Corrected=Frame_Information;
% Reset the strutures that will be combined
LocalizationsFinal={};
Frame_Information={};

for ijk=1:max(Came_from_image)
    LocalizationsFinal{ijk}=[];
    Frame_Information{ijk}=[];
end

for ksu=1:length(cut1array)
    ksu
    pp=.2;
    temp1=0;
    temp2=0;
    
    X1=Final_Loc_Blinking_Corrected{ksu}(:,1);
    X2=Final_Loc_Blinking_Corrected{ksu}(:,2);
    X3=Final_Loc_Blinking_Corrected{ksu}(:,3);
    X4=Final_Frame_Blinking_Corrected{ksu};
    
    
    addon=addonarray(ksu);
    cut1=cut1array{ksu};
    cut2=cut2array{ksu};
    cut3=cut3array{ksu};
    
    if cut3(2)==0
        cut3=[-1000,1000];
    end
    
     if ~isempty(Cutsstore)
    ssde=[cut1(1),cut2(1)]-Cutsstore;
    ssde2=[cut1(2),cut2(2)]-Cutsstore;
    if sum((sum(ssde')==0))==0 && sum((sum(ssde2')==0))==0
        IND=find(X1>=cut1(1) & X1<=cut1(2) & X2>=cut2(1) & X2<=cut2(2) & X3>=cut3(1) & X3<=cut3(2));
    elseif sum((sum(ssde')==0))==0
         IND=find(X1>=cut1(1) & X1<cut1(2) & X2>=cut2(1) & X2<cut2(2) & X3>=cut3(1) & X3<=cut3(2));
    elseif sum((sum(ssde2')==0))==0
         IND=find(X1>cut1(1) & X1<=cut1(2) & X2>cut2(1) & X2<=cut2(2) & X3>=cut3(1) & X3<=cut3(2));
    else
        IND=find(X1>cut1(1) & X1<cut1(2) & X2>cut2(1) & X2<cut2(2) & X3>=cut3(1) & X3<=cut3(2));
    end
    else
         IND=find(X1>=cut1(1) & X1<=cut1(2) & X2>=cut2(1) & X2<=cut2(2) & X3>=cut3(1) & X3<=cut3(2));
    end
        
        
    Cutsstore=[Cutsstore; [cut1(2), cut2(2)]];
    Cutsstore=[Cutsstore; [cut1(1), cut2(1)]];
    figure(1)
    scatter3(X1(IND),X2(IND),X3(IND),10,(X4(IND)),'filled')
    axis equal
    drawnow
    colormap jet
    pause(.5)

    
    LocalizationsFinal{Came_from_image(ksu)}=[LocalizationsFinal{Came_from_image(ksu)}; [X1(IND),X2(IND),X3(IND)]];
    app=X4(IND);
    Frame_Information{Came_from_image(ksu)}=[Frame_Information{Came_from_image(ksu)}; app(:)];
    
    
end

% Finally we go through and update the structures containing the data and
% save out the structures. With the needed information to reconstruct the
% images appropriatly.

Final_Loc_Blinking_Corrected=LocalizationsFinal;
Final_Frame_Blinking_Corrected=Frame_Information;
save(['Combined_Final_',Condition],'Final_Loc_Blinking_Corrected','Final_Frame_Blinking_Corrected')


%%
scatter(Final_Loc_Blinking_Corrected{1}(:,1),Final_Loc_Blinking_Corrected{1}(:,2),20, Frame_Information{1},'filled')
axis equal
    drawnow
    colormap jet