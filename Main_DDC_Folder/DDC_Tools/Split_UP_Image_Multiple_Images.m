%This is going to be a script to split up an image into multiple images so
%that DDC can act on them in //.

%The target number of localizations is 4000 (or whatever min_loc is)
%and there will be some buffer added to each side of the images to avoid boundary
%effects. -CHB 2019

clear

TrueLocalizations=[];
Photons={};
Photons{1}=[];
Resolution={};
min_loc=5000;

%Here you will be able to select the image you want to split up.
%This will allow you to pick out your own file, and will then load the data

[filename, pathname] = uigetfile({'*.mat'}, 'Select HMM .mat file');
if ( filename == 0 )
    disp('Error! No (or wrong) file selected!')
    return
end

%This will load in all of the localizations from your structure file
full_filename = [ pathname, filename ];
load(full_filename);

Condition=filename;


%If photons cell does not exist, generate one
%See user guide for specifics on this cell
if isempty(Photons{1})
    for i=1:length(Frame_Information)
        vv1=ones(1,length(Frame_Information{i}));
        Photons{i}=vv1(:);
    end
end


%These will temp store the split up images. 
LocalizationsFinal_Split={};
TrueLocalizations_Split={};
Frame_Information_Split={};

%I am adding the photons
Photons_Split={};

Came_from_image=[];
Parameters_to_split={};


%Just in case you do not actually know the true localizations, this way DDC
%can still work
if isempty(TrueLocalizations)
    TrueLocalizations={};
    for ijk=1:length(Frame_Information)
        TrueLocalizations{ijk}=[];
    end
end

%If the localizations are in 2D adjust the analysis so that the third
%colomn is just 0. This way DDC still works correctly as well as the
%methodology for splitting up the images.

for sdfv=1:length(LocalizationsFinal)
    if min(size(LocalizationsFinal{sdfv}))<3
        LocalizationsFinal{sdfv}(:,3)=LocalizationsFinal{sdfv}(:,2)*0;
        if min(size(TrueLocalizations{sdfv}))<3 && min(size(TrueLocalizations{sdfv}))>1
            TrueLocalizations{sdfv}(:,3)= TrueLocalizations{sdfv}(:,2)*0;
        end
    end
end

%Here we will go through an actually split up the images so that
%approximatly min_loc is within each image.
counter=0;
temp_numb_of_loc=[];
addonarray=[];
cut1array={};
cut2array={};
cut3array={};
%%
for ksu=1:length(Frame_Information)
    ksu
    pp=.2;
    temp1=0;
    temp2=0;
    
    X1=LocalizationsFinal{ksu}(:,1);
    X2=LocalizationsFinal{ksu}(:,2);
    X3=LocalizationsFinal{ksu}(:,3);
    
    if isempty(TrueLocalizations{ksu})==0
        X1t=TrueLocalizations{ksu}(:,1);
        X2t=TrueLocalizations{ksu}(:,2);
        X3t=TrueLocalizations{ksu}(:,3);
    end
    

    %This will make it so that if you only have 2d data you can still do the
    %split correctly
    if isempty(X3)
        X3=X2.*0;
    end
    
   
    X4=Frame_Information{ksu};
    
    addon=250; %Buffer region for the images when they get split up.
    
    scorefinal=Inf;
    onwers=1;
    while onwers<300
        
        disp(['Frac done with phase space search to split image' num2str(onwers/300)])

        temp_numb_of_loc=[];
        onwers=onwers+1;
        pp=rand*.95+.05;
        pp2=rand*.95+.05;
        
        
        %This will limit the phase space search if the cordinates are
        %actually in 2D. 
        if max(X3)~=0
            pp3=rand*.95+.05;
        else
            pp3=1;
        end
        
        
        
        
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
                    
                    IND=find(X1>=cut1(1)-addon & X1<=cut1(2)+addon & X2>=cut2(1)-addon & X2<=cut2(2)+addon & X3>=cut3(1)-addon & X3<=cut3(2)+addon);
                    
                    temp_numb_of_loc(end+1)=length(IND);
                end
            end
        end
        
        
        %This is the scoring function that we want to minimize, so that we
        %obtain the optimum number of localizations within our images. 
        Resid=(temp_numb_of_loc-min_loc);
        score=sum((Resid).^2)/length(temp_numb_of_loc);
        
        
        if score<scorefinal
            
            scorefinal=score;
            ppf=pp;
            ppf2=pp2;
            ppf3=pp3;
            onwers=0-onwers;
            temp_numb_of_locf=temp_numb_of_loc;
            
        end
        
    end
    
    
    
    

    addont=addon;
    
    %Now that we know the optimum way to split up our specific image, we
    %now split it up. 
    
    
    %%
%      ppf=.05;
%     ppf2=.05;
%     ppf3=1;
    
    pp=ppf;
    pp2=ppf2;
    pp3=ppf3;
    
    if length(X1)<min_loc
        
        ppf=1;
        ppf2=1;
        ppf3=1;
        
        pp=ppf;
        pp2=ppf2;
        pp3=ppf3;
        
    end
    
    flag=0;
    
    for i=1:ceil(1/pp)
        
        if flag==1
            break
        end
        
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
                
                addon=addont;
                
                IND=[];
                if length(X1)>min_loc
                    
                    while length(IND)<min_loc
                        addon= addon+10;
                        IND=find(X1>cut1(1)-addon & X1<cut1(2)+addon & X2>cut2(1)-addon & X2<cut2(2)+addon & X3>=cut3(1)-addon & X3<=cut3(2)+addon);
                    end
                   
                else
                    IND=1:length(X1);
                end
                
                
                while length(IND)>min_loc
                    addon= addon-10;
                    IND=find(X1>cut1(1)-addon & X1<cut1(2)+addon & X2>cut2(1)-addon & X2<cut2(2)+addon & X3>=cut3(1)-addon & X3<=cut3(2)+addon);
                  
                    if addon<150
                        break
                    end
                    
                end
                
                if isempty(TrueLocalizations{ksu})==0
                    INDt=find(X1t>cut1(1)-addon & X1t<cut1(2)+addon & X2t>cut2(1)-addon & X2t<cut2(2)+addon & X3t>=cut3(1)-addon & X3t<=cut3(2)+addon);
                end
                
                figure(1)
                scatter3(X1(IND),X2(IND),X3(IND),10,(X4(IND)),'filled')
                axis equal
                drawnow
                colormap jet
                pause(.5)
                counter=counter+1;
                
                LocalizationsFinal_Split{counter}=[X1(IND),X2(IND),X3(IND)];
                Photons_Split{counter}=Photons{ksu}(IND);
                
                if isempty(TrueLocalizations{ksu})==0
                    TrueLocalizations_Split{counter}=[X1t(INDt),X2t(INDt),X3t(INDt)];
                end
                
                cut1array{counter}=cut1;
                cut2array{counter}=cut2;
                cut3array{counter}=cut3;
                Frame_Information_Split{counter}=X4(IND);
                Came_from_image(counter)=ksu;
                Parameters_to_split{counter}=[ppf,ppf2,ppf3];
                addonarray(counter)=addon;
                temp_numb_of_loc(end+1)=length(IND);
                

                %This will account for images that allready have less than
                %the need number of localizations.
                if length(X1)<min_loc
                    flag=1;
                    break
                    
                end
            end
        end
    end
end

%%
%Finally we go through and update the structures containing the data and
%save out the structures. With the needed information to reconstruct the
%images appropriatly.

LocalizationsFinal=LocalizationsFinal_Split;
Frame_Information=Frame_Information_Split;
TrueLocalizations=TrueLocalizations_Split;
Photons=Photons_Split;

save(['Split_',Condition],'LocalizationsFinal','Frame_Information','addonarray','Parameters_to_split','Came_from_image','TrueLocalizations','cut1array','cut2array','cut3array','Resolution','Photons')

