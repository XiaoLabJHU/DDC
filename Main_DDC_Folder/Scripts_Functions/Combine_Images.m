%This will be a script to combine back the split images, so that the
%localizations can be observed within the cell like geometry. 

clear

%Just click run and click on the analyzed data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NOTE: this will only work on data that has been processed with DDC.

[filename, pathname] = uigetfile({'*.mat'}, 'Select HMM .mat file');
if ( filename == 0 )
    disp('Error! No (or wrong) file selected!')
    return
end

%This will load in all of the localizations from your structure file
full_filename = [ pathname, filename ];
load(full_filename);
load(filename(7:end))


Condition=filename;

pp=.2;
temp1=0;
temp2=0;
counter=1;
X1=LocalizationsFinal{1}(:,1);
X2=LocalizationsFinal{1}(:,2);
X3=LocalizationsFinal{1}(:,3);
X4=Frame_Information{1};


Combined_Image=[];
Combined_Image_Frame=[];


%Make sure you change this Chris
    pp=ppf
    pp2=ppf2
    pp3=ppf3
    
    
    counterss=0;
    for i=1:ceil(1/pp)
        if i*pp<1
            cut1=quantile(X1,[(i-1)*pp i*pp]);
        else
            cut1=quantile(X1,[(i-1)*pp 1]);
        end
        
        for ii=1:ceil(1/pp2)
            if ii*pp2<1
                cut2=quantile(X2,[(ii-1)*pp2 ii*pp2]);
            else
                cut2=quantile(X2,[(ii-1)*pp2 1]);
            end
            
            
            for iii=1:ceil(1/pp3)
                if iii*pp3<1
                    cut3=quantile(X3,[(iii-1)*pp3 iii*pp3]);
                else
                    cut3=quantile(X3,[(iii-1)*pp3 1]);
                end
                
                counterss=counterss+1;
        
                IND=find(Final_Loc_Blinking_Corrected{counterss}(:,1)>cut1(1) & Final_Loc_Blinking_Corrected{counterss}(:,1)<cut1(2) & Final_Loc_Blinking_Corrected{counterss}(:,2)>cut2(1) & Final_Loc_Blinking_Corrected{counterss}(:,2)<cut2(2) & Final_Loc_Blinking_Corrected{counterss}(:,3)>=cut3(1) & Final_Loc_Blinking_Corrected{counterss}(:,3)<=cut3(2));
                length(IND)
                Combined_Image=[Combined_Image; Final_Loc_Blinking_Corrected{counterss}(IND,1:end)];
                Combined_Image_Frame=[Combined_Image_Frame, Final_Frame_Blinking_Corrected{counterss}(IND)];

            end
        end
    end

 
       
save(['Combined_',Condition],'Combined_Image','Combined_Image_Frame')

%%
figure(1)

scatter(Combined_Image(:,1),Combined_Image(:,2),25,Combined_Image_Frame,'filled')
hold on 
plot(TrueLocalizations{1, 1}(:,1),TrueLocalizations{1, 1}(:,2),'r.')
