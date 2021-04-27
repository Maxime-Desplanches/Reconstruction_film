function [B] = vrai_carte(spos,nbs,type,A,kjq)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
B=zeros(size(A));
for i=1:nbs
    B(spos(i,1),spos(i,2))=10*normpdf(0,0,15);
end

switch type
    
    case 'normal'
        if nargin<5
            kjq=15;
        end
        for i=1:nbs
            B(spos(i,1),spos(i,2))=10*normpdf(0,0,kjq/5/1.5+10);
        end
        for k=1:nbs
            for i=1:size(B,1)
                for j=1:size(B,2)
                    
                    if ~isequal([i,j],spos(k,:))
                        
                        B(i,j)=B(i,j)+10*normpdf(sqrt((i-spos(k,1))^2+(j-spos(k,2))^2),0,kjq/5/1.5+10);
                        
                    end
                end
            end
        end
        
    case 'exp'
        if nargin<5
            kjq=15;
        end
        for i=1:nbs
            B(spos(i,1),spos(i,2))=10*exp(0);
        end
        for k=1:nbs
            for i=1:size(B,1)
                for j=1:size(B,2)
                    if ~isequal([i,j],spos(k,:))
                        
                        B(i,j)=B(i,j)+10*exp(-sqrt(( i - spos(k,1) )^2+(j-spos(k,2))^2)/4);
                        
                    end
                    
                end
            end
        end
    case 'diffuse'
        if nargin<5
            kjq=15;
        end
        B=imgaussfilt(B,10*kjq);
end

