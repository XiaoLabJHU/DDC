%% Combine_Images_Multiple_Images.m

%This script recombines your images after DDC has been run

% The target number of localizations is  min_loc and some buffer will be
% added to each side of the images to avoid boundary effects.
clear
files92307898=dir('An*');

Cutsstore=[];
for iwen=1:length(files92307898)
    
filename5578599999=files92307898(iwen).name;

load(filename5578599999);

Condition=filename5578599999;

%Final_Loc_Blinking_Corrected=LocalizationsFinal;
%Final_Frame_Blinking_Corrected=Frame_Information;
%Reset the strutures that will be combined


%For making the thresholded images
% for isnese=1:length(LocalizationsFinal)
% loc=LocalizationsFinal{isnese};
% fram=Frame_Information{isnese};
% [LocalizationsFinal{isnese}, Frame_Information{isnese}]=Threshold(loc, fram, 40, 100);
% end
% 
 if length(Final_Loc_Blinking_Corrected)~=0
% Final_Loc_Blinking_Corrected=LocalizationsFinal;
%Final_Frame_Blinking_Corrected=Frame_Information;

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
%     figure(1)
%     scatter3(X1(IND),X2(IND),X3(IND),10,(X4(IND)),'filled')
%     axis equal
%     drawnow
%     colormap jet
%     pause(.5)

    
    LocalizationsFinal{Came_from_image(ksu)}=[LocalizationsFinal{Came_from_image(ksu)}; [X1(IND),X2(IND),X3(IND)]+[randi([-2,2],length(IND),1),randi([-2,2],length(IND),1),randi([-2,2],length(IND),1)]];
    app=X4(IND);
    Frame_Information{Came_from_image(ksu)}=[Frame_Information{Came_from_image(ksu)}; app(:)];
    
    
end

% Finally we go through and update the structures containing the data and
% save out the structures. With the needed information to reconstruct the
% images appropriatly.

Final_Loc_Blinking_Corrected=LocalizationsFinal;
Final_Frame_Blinking_Corrected=Frame_Information;
save(['DDC_Combined_Final_',Condition],'Final_Loc_Blinking_Corrected','Final_Frame_Blinking_Corrected')
%save(['DDC_Combined_Final_',Condition],'Final_Loc_Blinking_Corrected','Final_Frame_Blinking_Corrected','cut5','cut6')

%%
% Condition
% %length(Final_Loc_Blinking_Corrected{1}(:,1))/length(LocalizationsFinal
% figure(1)
% scatter(Final_Loc_Blinking_Corrected{1}(:,1),Final_Loc_Blinking_Corrected{1}(:,2),20, Frame_Information{1},'filled')
% axis equal
%     drawnow
%     colormap jet
%   
%     
%loc1=LocalizationsFinal{1};
%loc=Final_Loc_Blinking_Corrected{1};
%fram=Frame_Information{v};
%[loc, Thresh_frame]=Threshold(loc, fram, Thresh_Time, Thresh_dist);
%[denstiy_map1]=PALMplot(loc(:,1:2), loc1(:,1:2));

%figure
%imshow(((denstiy_map1)))
%title(Condition)
%caxis([0 .2])
%colormap hot
%drawnow
 end
end




%%
figure
loc1=Final_Loc_Blinking_Corrected{1};
loc=Final_Loc_Blinking_Corrected{1};
%loc=LocalizationsFinal{1};
%fram=Frame_Information{v};
%[loc, Thresh_frame]=Threshold(loc, fram, Thresh_Time, Thresh_dist);
[denstiy_map1]=PALMplot(loc(:,1:2), loc1(:,1:2));

figure
imagesc(((denstiy_map1)))
colormap jet
caxis([0 2.5])
%%
%%


%%
figure
histogram(denstiy_map1(denstiy_map1>0),5000)

%%
%%
figure
%loc=LocalizationsFinal{1};
%fram=Frame_Information{1};
%loc1=Final_Loc_Blinking_Corrected{1};
loc=Final_Loc_Blinking_Corrected{1};
fram=Final_Frame_Blinking_Corrected{1};
%loc1=Final_Loc_Blinking_Corrected{1};
%loc=Final_Loc_Blinking_Corrected{1};
%fram=Final_Frame_Blinking_Corrected{1};
scatter(loc(:,1),loc(:,2),10,fram(:),'filled')
colormap jet
axis equal 

%%
loc=LocalizationsFinal{1};
fram=Frame_Information{1};
hold on 
plot(loc(:,1),loc(:,2),'+')