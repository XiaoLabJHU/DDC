%This is going to be a script that evaluates each of the images and
%compairs them to before and after the MCMC approach:) 

%First step: make sure you load in the data that you analyzed with DDC
%after you have loaded in your data just click run and it will generate the
%last histogram of the USER GUIDE.

%% 

stores=[];
stores2=[];
biner=0:.1:1;
for CEL=1:length(LocalizationsFinal)
    [mm,mmm]=max(Lik{1, CEL})
   stores(end+1)= RelScore{1, CEL}(1)
  % stores(end+1)=( abs(Numb_of_Loc{1, CEL}(mmm)-length(TrueLocalizations{CEL})))
  stores2(end+1)=RelScore{1, CEL}(mmm)
% stores2(end+1)=( abs(Numb_of_Loc{1, CEL}(1)-length(TrueLocalizations{CEL})))
end

%% Plotting Image Error

% See Figure 5 in USER GUIDE

figure(1)
histogram(stores,biner,'Normalization','prob')
mean(stores)
hold on 
histogram(stores2,biner,'Normalization','prob')
mean(stores2)
legend('Image Error No MCMC','Image Error After MCMC')
ylabel('Prob','FontSize',20)
xlabel('Image Error','FontSize',20)
title('Image Error')
set(gca,'FontSize',20)