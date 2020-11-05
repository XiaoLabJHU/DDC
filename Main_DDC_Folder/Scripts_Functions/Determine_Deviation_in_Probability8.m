function [Deviation_in_Probability, Prob_Distributions, Dscale_store]=Determine_Deviation_in_Probability8(True_Distribuiton, Distribution_for_Blink, bins, A, LocalizationsFinal, Frame_Information, loc, fram, Resolution,  True_Distribuiton2, Distribution_for_Blink2, bins2, X_overall, M_mat)

%X_overall=[];
%This is the script that makes the  probability of a localization being
%a blink of a SINGLE previous localization, a certain distance away with a particular frame difference.

%It also defines the distribution that are used in the calculation of the
%likilihood.

%This is the script that does the fitting!!!!

Prob_Distributions=[];

D = (pdist(LocalizationsFinal));

Z2 = pdist([[1:length(Frame_Information)]'*0,Frame_Information(:)]);


%The following commented out section was used to verify some of the
%calculations within DDC
%{
if length(loc)>1
    
    Df = pdist(loc);
    
    Z22 = pdist([[1:length(fram)]'*0 ,fram(:)]);
    
    [C ,ia,ib] = intersect(D,Df);
    
    D(ia)=[];
    
    Z2(ia)=[];
    
end
%}


Deviation_in_Probability=[];

Dscale_store=[];

D_Scale_Max=0;


if length(X_overall)==0
    
    for w=1:A
        
        D_Blink=D(Z2==w);
        
        Temp_Distribution=((histcounts(D_Blink,bins2,'Normalization','prob')));
        
        y=Temp_Distribution;
        
        t=[True_Distribuiton2; Distribution_for_Blink2];
        
        F=@(x,xdata) x.*xdata(1,:)+(1-x).*xdata(2,:);
        
        x0=1;
        opts = optimset('Display','off');
        [x,~,~,~,~] = lsqcurvefit(F,x0,t,y,[],[],opts);
        
        D_Scale=x;
        
        Dscale_store(end+1)=D_Scale;
        
    end
    
else
    Dscale_store=X_overall;
end


if isempty(loc)
    Dscale_store5=Dscale_store;
end




Dscale_store(Dscale_store>1)=1;

Dscale_store(Dscale_store<0)=.0000001;

Dscale_store(end)=1;

Deviation_in_Probability=[];




if ~isempty(M_mat)
    Deviation_in_Probability=M_mat;
else
    
    for w=1:A
        
        D_Blink=D(Z2==w);
        
        %  Temp_Distribution=((histcounts(D_Blink,bins,'Normalization','prob')));
        
        D_Scale=Dscale_store(w);
        
        Temp_Distribution2=Distribution_for_Blink*(1-D_Scale)+True_Distribuiton*(D_Scale);
        
        %See the equation in determining the probability section of the paper.
        Combined=(Temp_Distribution2-D_Scale*True_Distribuiton)./(Temp_Distribution2);
        
        Deviation_in_Probability(w,:)=Combined;
        
    end
    
    
end

disp('Done with First Step of Calculating Joint Probability')

%D_Scale2=D_Scale;

%Time to determine the number of localizations approx
Traj_num=randperm(length(LocalizationsFinal),length(LocalizationsFinal));

[loc, ~]=Eliminate_Blinking_De_Loc15_MCMC(LocalizationsFinal, Frame_Information, Resolution, A, Deviation_in_Probability, Traj_num);

nb=(length(LocalizationsFinal)-sum(loc))/sum(loc);

D_Scale_Max=0;

for w=1:A
    
    D_Blink=D(Z2==w);
    
    Temp_Distribution=((histcounts(D_Blink,bins2,'Normalization','prob')));
    
    y=Temp_Distribution;
    
    %I have modified this, so that it is more correct!!!!
    
    y=y-Dscale_store(w)*(1-((nb+nb+nb*nb)/(1+(nb+nb+nb*nb)))).*True_Distribuiton2;
    
    y(y<0)=0;
    
    y=y/(sum(y));
    
    t=[True_Distribuiton2; Distribution_for_Blink2];
    
    F=@(x,xdata) x.*xdata(1,:)+(1-x).*xdata(2,:);
    
    x0=1;
    opts = optimset('Display','off');
    [x,~,~,~,~] = lsqcurvefit(F,x0,t,y,[],[],opts);
    
    D_Scale=x;
    
    Dscale_store(w)=D_Scale;
    
end


Dscale_store(Dscale_store>1)=1;
Dscale_store(Dscale_store<0)=.0000001;
Dscale_store(end)=1;


for w=1:A
    
    D_Blink=D(Z2==w);
    
    % Temp_Distribution=((histcounts(D_Blink,bins,'Normalization','prob')));
    
    D_Scale=Dscale_store(w);
    
    Temp_Distribution2=Distribution_for_Blink*(1-D_Scale)+True_Distribuiton*(D_Scale);
    
    Prob_Distributions(w,:)=Temp_Distribution2;
    
    Combined=(Temp_Distribution2-D_Scale*True_Distribuiton)./(Temp_Distribution2);
    
end




disp('Done with First Step of Calculating Joint Probability')




%This is just in case one of the bins is zero, set it to be ten times lower
%than the lowest bin.

Prob_Distributions=Prob_Distributions+min(Prob_Distributions(Prob_Distributions>0))/10;




