function [cval] = creation_valeur_3d(cpos,spos,nbc,nbs,type,k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
cval=zeros(nbc,1);
switch type
    case 'normal'
        for i=1:nbs
            for j=1:nbc
                cval(j)=cval(j)+10*normpdf(sqrt(( cpos(j,1)-spos(i,1) )^2+(cpos(j,2)-spos(i,2))^2),0,k/5+10);
            end
        end
        
    case 'exp'
         for i=1:nbs
            for j=1:nbc
                cval(j)=cval(j)+10*exp(-sqrt(( cpos(j,1) - spos(i,1) )^2+(cpos(j,2)-spos(i,2))^2)/51-k);
            end
         end
         
    case 'diffuse'
        for i=1:nbs
            for j=1:nbc
                cval(j)=cval(j)+10*normpdf(sqrt(( cpos(j,1) - spos(i,1) )^2+(cpos(j,2)-spos(i,2))^2),0,20+k);
            end
        end  
    case 'erreur'
                for i=1:nbs
            for j=1:nbc
                cval(j)=cval(j)+10*normpdf(sqrt(( cpos(j,1)-spos(i,1) )^2+(cpos(j,2)-spos(i,2))^2),0,floor(k/3)+5)*(1+abs(randn(1))/20);
            end
        end
end


