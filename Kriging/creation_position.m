function [cpos,posl,posc] = creation_position(type,nbc,A,data)
if nargin==3
    data=[];
end
cpos=zeros(nbc,2);


switch type
    case 'random'
        nargout=1;
        i=1;
        while cpos(nbc,1)==0
            val=1+floor([size(A,1) size(A,2)].*rand(1,2));
            if isempty(find(ismember(cpos,val,'rows')))
                cpos(i,1)=val(1);
                cpos(i,2)=val(2);
                i=i+1;
            end
        end
    case 'source_centre'
        nargout=1;
        i=1;
        while cpos(nbc,1)==0
            val=1+floor([size(A,1)/2 size(A,2)/2].*rand(1,2));
            if isempty(find(ismember(cpos,val,'rows')))
                cpos(i,1)=val(1)+floor(size(A,1)/4);
                cpos(i,2)=val(2)+floor(size(A,2)/4);
                i=i+1;
            end
        end
    case 'linear'
        nargout=1;
        ecart=size(A,1)*size(A,2)/nbc;
        for i=1:nbc
            cpos(i,1)=floor(mod(i*ecart,size(A,1)))+1;
            cpos(i,2)=floor(i*ecart/size(A,1))+1;
        end
        if abs(mean(cpos(:,1))-size(A,1)/2)>1
            a=floor(abs(max(cpos(:,1))-size(A,1))/2);
            for i=1:nbc
                cpos(i,1)=cpos(i,1)+a;
            end
        end
    case 'connu'
        nargout=1;
        X=load(data);
        cpos=X.ind_sources;
        
    case 'grille'
        [j,k]=size(A);
        X=true;
        m=floor(sqrt(nbc));
        n=floor(sqrt(nbc))+1;
        while X
            if mod(nbc,m)==0
                X=false;
            else
                m=m-1;
            end
        end
        if j>k
            n=nbc/m;
        else
            n=m;
            m=nbc/n;
        end
        posl=rand(1,n);
        posc=rand(1,m);
        vall=floor(j/n);
        valc=floor(k/m);
        for i=1:max(n,m)
            if i<=m
                posc(i)=floor((posc(i)+i-1)*valc)+1;
            end
            if i<=n
                posl(i)=floor((posl(i)+i-1)*vall)+1;
            end
        end
        [O,P]=meshgrid(posl,posc);
        cpos=[O(:),P(:)];
    case 'grille_reg'
        [j,k]=size(A);
        X=true;
        m=floor(sqrt(nbc));
        n=floor(sqrt(nbc))+1;
        while X
            if mod(nbc,m)==0
                X=false;
            else
                m=m-1;
            end
        end
        if j>k
            n=nbc/m;
        else
            n=m;
            m=nbc/n;
        end
        posl=zeros(1,n)+0.5;
        posc=zeros(1,m)+0.5;
        vall=floor(j/n);
        valc=floor(k/m);
        for i=1:max(n,m)
            if i<=m
                posc(i)=floor((posc(i)+i-1)*valc)+1;
            end
            if i<=n
                posl(i)=floor((posl(i)+i-1)*vall)+1;
            end
        end
        [O,P]=meshgrid(posl,posc);
        cpos=[O(:),P(:)];
    case 'source_proche'
        nargout=1;
        
        val=1+floor([size(A,1)/2 size(A,2)/2].*rand(1,2));
        cpos(1,1)=val(1)+floor(size(A,1)/4);
        cpos(1,2)=val(2)+floor(size(A,2)/4);
        
        for i=2:size(cpos,1)
            cpos(i,1)=cpos(i-1,1)+5;
            cpos(i,2)=cpos(i-1,2)+5;
        end
    case 'source_eloigne'
        [I,J]=size(A);
        nargout=1;
        for i=1:nbc
        cpos(i,1)=floor(i*(I-I/4)/nbc+I/8);
        cpos(i,2)=floor(i*(J-J/4)/nbc+J/8);

        end
end

