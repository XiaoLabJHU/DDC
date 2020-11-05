%THis is going to be a script to combine multiple images
clear
files=dir('An*');

for apple=1:length(files)
    files=dir('An*');
    ttse=files(apple).name;
    load(files(apple).name)
    
    
    Condition=ttse;

% Final_Loc_Blinking_Corrected=LocalizationsFinal;
% Final_Frame_Blinking_Corrected=Frame_Information;
%Here we are going to reset the strutures as we are going to combine
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
    
    IND=find(X1>=cut1(1) & X1<=cut1(2) & X2>=cut2(1) & X2<=cut2(2) & X3>=cut3(1) & X3<=cut3(2));
    
%     figure(1)
%     scatter3(X1(IND),X2(IND),X3(IND),10,(X4(IND)),'filled')
%     axis equal
%     drawnow
%     colormap jet
%     pause(.5)

    
    LocalizationsFinal{Came_from_image(ksu)}=[LocalizationsFinal{Came_from_image(ksu)}; [X1(IND),X2(IND),X3(IND)]];
    app=X4(IND);
    Frame_Information{Came_from_image(ksu)}=[Frame_Information{Came_from_image(ksu)}; app(:)];
    
    
end

%%
%Finally we go through and update the structures containing the data and
%save out the structures. With the needed information to reconstruct the
%images appropriatly.

Final_Loc_Blinking_Corrected=LocalizationsFinal;
Final_Frame_Blinking_Corrected=Frame_Information;
clf 
indss=randperm(length(Final_Loc_Blinking_Corrected{1}(:,1)));
scatter(Final_Loc_Blinking_Corrected{1}(indss,1),Final_Loc_Blinking_Corrected{1}(indss,2),10,Final_Frame_Blinking_Corrected{1}(indss),'filled')
colormap jet
axis equal

pause
save(['Combined_Final_',Condition],'Final_Loc_Blinking_Corrected','Final_Frame_Blinking_Corrected')
end
%%