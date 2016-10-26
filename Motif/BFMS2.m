% function 2 --BFMS2
function [cstr sc pos] = BFMS2(DNA,lmer,iftrace)
	disp('Algorithm 2 - Branch and Bound Motif Search');
	[t n] = size(DNA);
	s=ones(1,t);
	bestscore=0;
	i=1;
	while i>0
		if i<t
            %optimisticScore=score(s,DNA,i)+(t-i)*lmer;
			optimisticScore=score(s,DNA,lmer)+(t-i)*lmer;
			if optimisticScore<bestscore
				if iftrace == 1
					sol = '';
					for k=1:i
						sol = strcat(sol, sprintf(' %d', s(k)));
					end
					for k=i+1:t
						sol = strcat(sol, ' -');
					end
					disp(sprintf('Bypass candidate: (%s )', sol));
				end
				[s,i]=Bypass(s,i,t,n-lmer+1);
			else
				[s,i]=nextvertex(s,i,t,n-lmer+1);
			end
        else
            %if score(s,DNA,t)>bestsocre
            if score(s,DNA,lmer)>bestscore
                [sc cstr]=score(s,DNA,lmer);
				%[sc cstr]=score(s,DNA,t);
				bestscore =sc;
				bestmotif=cstr;
                pos=s;
                sol='';
            	for i =1:length(s)
                    sol = strcat(sol,sprintf(' %d',s(i)));
            	end
            disp(sprintf('\tUpdate: bestscore=%3d  at (%s )', bestscore,sol));
            end
            [s,i]=nextvertex(s,i,t,n-lmer+1);
        end
    end
end
function [out level]=Bypass(a,i,L,k) %a:array of digits i:prefix length L;maximum length k:max digit value
	for j=i:-1:1
		if a(j)<k
			a(j)=a(j)+1;
			level=j;
			out =a(1:level); 
			return
		end
	end
	out=a;
	level=0;
	return
end
function [out level]=nextvertex(a,i,L,k)
    if i<L
        a(i+1)=1;
        out =a;
        level=i+1;
        return 
    else
        for j=L:-1:1
            if a(j)<k
                a(j)=a(j)+1;
            	level=j;
                out =a(1:level);
                return 
            end
        end
    end
    out=a;
    level=0;
    return
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



