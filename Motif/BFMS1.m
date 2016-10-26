% function 1 --BFMS1
function [cstr sc pos] = BFMS1(DNA,lmer,iftrace)
    disp('Algorithm 1 - Simple Brute Forcef Motif Search');
    [t n] = size(DNA);
    s =ones(1,t);
    bestscore = score(s,DNA,lmer);
    sol = '';
    for i =1:length(s)
        sol = strcat(sol,sprintf(' %d',s(i)));
    end
    disp(sprintf('\tUpdate: bestscore=%3d  at (%s )', bestscore,sol));
    while true
        s=nextleaf(s,t,n-lmer+1); %n-lmer+1 is the maxinum number for si
        newscore = score(s,DNA,lmer);
        [sc cstr]=score(s,DNA,lmer);
        if iftrace == 1
            if s(t-1) ==1 && s(t) ==1
                sol = '';
                for i =1:length(s)
                    sol = strcat(sol,sprintf(' %d',s(i)));
                end
                disp(sprintf('Pass candidate: (%s )', sol));
            end
        end
        if newscore > bestscore
            bestscore=newscore;
            bestMotif=cstr;
            sol='';
            for i =1:length(s)
                    sol = strcat(sol,sprintf(' %d',s(i)));
            end
            disp(sprintf('\tUpdate: bestscore=%3d  at (%s )', bestscore,sol));
        end
         %if s=one(lmer)
        pos=s;
           % [sc cstr]=score(s,DNA,lmer);
        % end
    end
end
function [sc cstr]=score(s,DNA,l)
    prof=profile_dna(s,DNA,l);
    [sc cstr]=score_p(prof);
end

%profile_dna.m file
function p = profile_dna(s,DNA,l)
    for i =1:length(s)
         align(i,1:l)=DNA(i,s(i):s(i)+l-1);
         %p(i,1:l)=atgc2num(align(i,1:l));
    end
    for j=1:l
        pj=align(:,j);
        p(1,j)=length(find(pj=='a'|pj=='A'));
        p(2,j)=length(find(pj=='c'|pj=='C'));
        p(3,j)=length(find(pj=='g'|pj=='G'));
        p(4,j)=length(find(pj=='t'|pj=='T'));
    end       
end

%score_p.m
function [sc cstr]=score_p(profile)
    [row,col]=size(profile);
    for i=1:col
        p=find(profile(:,i)==max(profile(:,i)));
        cstr1(i)=p(1);
        sc1(i)=max(profile(:,i));
    end
    sc =0;
    for j =1:col
        sc=sc+sc1(j);
    end
    cstr=num2atgc(cstr1);
end
function str = num2atgc(mer)
    for i =1:length(mer)
        if mer(i) == 1
            str(i)='a';
        elseif mer(i) ==2
            str(i)='t';
        elseif mer(i) ==3 
            str(i)='g';
        elseif mer(i) ==4 
            str(i)='c';
        else
            return
        end 
    end
end
function leaf = nextleaf(a,L,k)
    for i =L:-1:1
        if a(i)<k
            a(i) = a(i) +1;
            leaf =a;
            return ;
        else
            a(i) =1;
        end
    end
    leaf =a;
end

            