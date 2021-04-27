function [cval] = creation_valeur(cpos,spos,nbc,nbs,type)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
cval=zeros(nbc,1);
switch type
    case 'normal'
        for i=1:nbs
            for j=1:nbc
                cval(j)=cval(j)+10*normpdf(sqrt(( cpos(j,1)-spos(i,1) )^2+(cpos(j,2)-spos(i,2))^2),0,15);
            end
        end
        
    case 'exp'
         for i=1:nbs
            for j=1:nbc
                cval(j)=cval(j)+10*exp(-sqrt(( cpos(j,1) - spos(i,1) )^2+(cpos(j,2)-spos(i,2))^2)/4);
            end
         end
         
    case 'diffuse'
        for i=1:nbs
            for j=1:nbc
                cval(j)=cval(j)+10*normpdf(sqrt(( cpos(j,1) - spos(i,1) )^2+(cpos(j,2)-spos(i,2))^2),0,100);
            end
        end         
end

