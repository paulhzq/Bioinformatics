% PDP problem


function PDP(L, iftrace)
	
	L = sort(L);	% make L sorted for easy implementation
	n = length(L); % n is the length of L

	width = max(L); % width is the maximum element of L
	L(L == width) = [ ];	% delete(width, L)
	X = [0 width]; 

	PLACE(L, X, width, 0, iftrace);
end

    
function PLACE(L, X, width, level, iftrace)

    indent = '----';
    for i=1:level
        indent = strcat(indent, '----');
    end
    
    if(isempty(L))% if L is empty.
		sol = '';
		XX = sort(X);
        for i=1:length(X)
			sol = strcat(sol, sprintf(' %d', XX(i)));
        end
		disp(sprintf('%s Solution found:%s', indent, sol));
    else
        y = max(L);   % y is the maximum element in L
		dy = abs(y-X);  % dy is delta(y, X)
		z = abs(width-y);
		dz = abs(z-X);  % dz is delta(width-y, X)
        
        if allin(dy, L) % delta(y,X) in L
           if iftrace == 1
               disp(sprintf('%s try   y=%d',indent,y));
           end
           X=[X,y];% add y to X
           L=Remove(dy,L);% remove all delta(y,X) in L
           PLACE(L,X,width, level+1, iftrace);
           X= Remove(y,X); %remove y from X
           L = sort([dy L]); % add delta(y,X) to L
        else
            if iftrace == 1
                disp(sprintf('%s try y=%d FAILS', indent,y));
            end
        end
        
        if allin(dz, L) % delta(width-y,X) in L,it's pretty much the same as allin(dy,L)
            if iftrace == 1
				disp(sprintf('%s try w-y=%d', indent, z));
            end
			X = [X,z] ;	  % add width-y to X
			L = Remove(dz,L); % remove all delta(width-y) in L
			PLACE(L, X, width, level+1, iftrace);
			X = Remove(z,X) ;	% remove width-y from X
			
			L = sort([L,dz]);	% add delta(width-y,X) to L
		else
			if (iftrace == 1)
				disp(sprintf('%s try w-y=%d FAILS', indent, z));
            end
  
        end
        
    end
end

    
        
        
function OUT = Remove(input_arr1,input_arr2)
    tmp = input_arr2;
    for n = 1:length(input_arr1)
        [is_not,col] = ismember(input_arr1(n),tmp);
        if(is_not == 1)
            tmp(:,col) = [];
        end
    end
    OUT = tmp;
end

function bool = allin(input_arr1,input_arr2)
    tmp = input_arr2;
    for n = 1:length(input_arr1)
        [is_not,col] = ismember(input_arr1(n),tmp);
        if(is_not == 0)
            bool = 0;
            break
        else if (is_not == 1)
                tmp(:,col) = [];
            end
            bool = 1;
        end
    end
end

 
        
            
        
    
    
    
    