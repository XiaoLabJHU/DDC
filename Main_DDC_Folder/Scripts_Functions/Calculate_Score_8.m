   

% 
% %Now here we will calculate our score
% lik=0;
%     for i=1:A
%         if length(DistanceMatrix3{i}(:))>1
%            if length(DistancMatrix3{i})<numb_of_dat
%             lik =lik+ log(pdf(pd,mean(DistanceMatrix3{i}(:))));
%            else
%             lik =lik+log(pdf(pd, mean(DistanceMatrix3{i}(randperm(length(DistancMatrix3{i}),numb_of_dat)))));
%            end
%         end
%     end
%     
%     
    

D = pdist(loc);
Z22 = pdist([[1:length(fram)]'*0 ,fram(:)]);
D=D(Z22<A);
Z22=Z22(Z22<A);
D=D(Z22>0);


D_temp=D_Overall;
Z2_temp=Z2;
%Z2_temp(Z2_temp>100)=100;


[C ,ia,ib] = intersect(D_Overall,D);
D_temp(ia)=[];
Z2_temp(ia)=[];
D = floor(D/Resolution);
D_temp=floor(D_temp/Resolution);
%length(D)+length(D_temp);
D(D+1>size(Prob_Distributions,2))=size(Prob_Distributions,2)-1;
D_temp(D_temp+1>size(Prob_Distributions,2))=size(Prob_Distributions,2)-1;

lik =sum(Prob_Distributions(end,D+1));
lik2=0;
for iin=1:A
    if ~isempty(D_temp(Z2_temp==iin))
        lik2=lik2+sum(Prob_Distributions(iin,D_temp(Z2_temp==iin)+1));
%     else
%         if length(D_temp(Z2_temp==iin))==0
%         hey=1
%         end
    end
end

%290974
% D3=[];
% for oi=1:10
% D2 = round(DistanceControl2(randperm(length(DistanceControl2),length(D))),-1);
% D2(D2==0)=10;
% D3(end+1) =mean(sum(Probabilities(D2/10+1)));  
% end
% lik2=mean(D3);
    
%   %Now here we will calculate our score
% lik=0;
% lik2=0;
% 
% distfinal=[];
%     for i=1:A
%         if length(DistanceMatrix3{i}(:))>1
%             distfinal=[DistanceMatrix3{i}];
%          
%             lik =lik+ log(pdf(pd_samples{length(distfinal)},[mean(distfinal),std(distfinal),moment(distfinal,3)]));
%            % lik2 =lik2+ log(pdf(pd_samples2{length(distfinal)},std(distfinal)));
%            
%         end
%     end