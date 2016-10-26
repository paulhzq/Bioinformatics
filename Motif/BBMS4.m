function [word,distance,post]=BBMS4(DNA,lmer,iftrace)
	disp('Algorithm 4 - Branch and Bound Median String Search');
	[t n] = size(DNA);
    s =ones(1,lmer); % s is the num to stand for the acgc.
	bestDistance=t*lmer;
	i=1;
	while i>0
		if i<lmer
            prefix =num2atgc(s);
            optimisticDistance=TotalDistance(prefix,DNA);
            if optimisticDistance>bestDistance
                if iftrace == 1
					sol = '';
					for i =1:length(prefix)
                        sol = strcat(sol,sprintf(' %s',prefix(i)));
                    end
                    for k=length(prefix)+1:lmer
						sol = strcat(sol, ' -');
					end
                    %disp(sprintf('Pass candidate: (%s )', sol));
					disp(sprintf('Bypass candidate: (%s )', sol));
				end
                
                [s i]=Bypass(s,i,lmer,4);
            else
                [s,i]=nextvertex(s,i,lmer,4);
            end
        else
			word_ini =num2atgc(s); % word ← nucleotide string corresponding to (s1, s2, . . .sl)
% 			if iftrace == 1
%                 if word_ini(4) =='C' && word_ini(3) =='C'&& word_ini(2)=='C'
%                     sol = '';
%                     for i =1:length(word_ini)
%                         sol = strcat(sol,sprintf(' %s',word_ini(i)));
%                     end
%                     disp(sprintf('Pass candidate: (%s )', sol));
%                 end
%             end
            if TotalDistance(word_ini,DNA)<bestDistance   %[dist,post]=TotalDistance(word,DNA)
				bestDistance=TotalDistance(word_ini,DNA);
				bestword = word_ini;
				word=bestword;
				[distance,post]=TotalDistance(word_ini,DNA);
				%pos=post;
				sol='';
            	for i =1:length(post)
                    sol = strcat(sol,sprintf(' %d',post(i)));
            	end
            disp(sprintf('\tUpdate: bestDist=%3d,bestWord=%s at (%s )', bestDistance,bestword,sol));
        	end
        	[s,i]=nextvertex(s,i,lmer,4);
        end
    end
end

               
%a=1 t=2,g=3,c=4
function mer = atgc2num(str)
    for i =1:length(str)
        if str(i) == 'a' | str(i) =='A'
            mer(i)=1;
        elseif str(i) =='t' | str(i) == 'T'
            mer(i)=2;
        elseif str(i) =='g' | str(i) =='G'
            mer(i)=3;
        elseif str(i) =='c' | str(i) =='C'
            mer(i)=4;
        else
            return
        end 
    end
end

function str = num2atgc(mer)
    for i =1:length(mer)
        if mer(i) == 1
            str(i)='A';
        elseif mer(i) ==2
            str(i)='T';
        elseif mer(i) ==3 
            str(i)='G';
        elseif mer(i) ==4 
            str(i)='C';
        else
            return
        end 
    end
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
function [dist,post]=TotalDistance(word,DNA)
	[t n]=size(DNA);
	for i=1:t
		[dist1(i) post(i)]=hm(word,DNA(i,:));
	end
	dist=0;
	for j=1:t
		dist=dist+dist1(j);
	end
end

function [min_dis pos]=hm(v,DNA_one_row)
	len = length(DNA_one_row);
	min = length(v);
	for i=1:len-length(v)+1
		outstr=substr(DNA_one_row,i,length(v));
		newdis=hamming_distance(v,outstr);
		if newdis<min
			min = newdis;
			min_dis=min;
			pos = i;
		end
	end
end

function dis = hamming_distance(v,d)
	len=length(v);
	total=0;
	for i=1:len
		if v(i)~=d(i)
			total=total+1;
		end
	end
	dis = total;
end

function outstr=substr(str,offset,len)
	outstr=str(offset:offset+len-1);
end

