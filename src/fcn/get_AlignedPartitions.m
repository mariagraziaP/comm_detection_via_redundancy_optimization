function new_vector = get_AlignedPartitions(C1,C2)

clust_number1 = length(unique(C1));
clust_number2 = length(unique(C2));

if clust_number1 >= clust_number2
    
    L1 = C1;
    L2 = C2;
    
    for i=1:clust_number1
        a = find(C1);
        num1(i) = C1(a(1));
        ind = find(C1==num1(i));
        clusterSize_C1(i) = length(ind);
        C1(ind) = 0;
    end
    clear a C1 ind i
    
    new_vector = zeros(size(L1,1),1);
    
    for i=1:clust_number2
        a = find(C2);
        num2(i) = C2(a(1));
        ind = find(C2==num2(i));
        clusterSize_C2(i) = length(ind);
        C2(ind) = 0;
    end
    clear a C2 ind i
    
    new_vector = zeros(size(L2));
    
    for k=1:clust_number1
        
        p = zeros(clust_number2,1);
        y = zeros(clusterSize_C1(k),1);
        clust = find(L1==num1(k));
        for j=1:clust_number2
            y(1:clusterSize_C1(k)) = num2(j);
            p(j) = sum(eq(y,L2(clust)));
        end
        
        maxcorr(k) = max(p);
        maxpos = find(p==max(p));
        maxcl(k) = num2(maxpos(1));
        
        posiz{k} = maxpos;
        freq(k,:) = p;
        
        clear p y
    end
    
    f1 = 1:clust_number1;
    f2 = 1:clust_number2;
    
    for w=1:clust_number2
        maxx = max(max(freq(f1,f2)));
        indmax = find(freq(f1,f2)==maxx);
        indmax = indmax(1);
        clL2 = ceil(indmax/length(f1));
        clL1 = indmax-length(f1)*floor(indmax/length(f1));
        if clL1==0
            clL1 = length(f1);
        end
        indL2 = find(L2==num2(f2(clL2)));
        new_vector(indL2) = num1(f1(clL1));
        f1(clL1) = [];
        f2(clL2) = [];
        
        clear indL2
    end
    
    
else
    L1 = C1;
    L2 = C2;
    
    for i=1:clust_number1
        a = find(C1);
        num1(i) = C1(a(1));
        ind = find(C1==num1(i));
        clusterSize_C1(i) = length(ind);
        C1(ind) = 0;
    end
    clear a C1 ind i
    
    new_vector = zeros(size(L1,1),1);
    
    for i=1:clust_number2
        a = find(C2);
        num2(i) = C2(a(1));
        ind = find(C2==num2(i));
        clusterSize_C2(i) = length(ind);
        C2(ind) = 0;
    end
    clear a C2 ind i
    
    new_vector = zeros(size(L2));
    
    for k=1:clust_number1
        
        p = zeros(clust_number2,1);
        y = zeros(clusterSize_C1(k),1);
        clust = find(L1==num1(k));
        for j=1:clust_number2
            y(1:clusterSize_C1(k)) = num2(j);
            p(j) = sum(eq(y,L2(clust)));
        end
        
        maxcorr(k) = max(p);
        maxpos = find(p==max(p));
        maxcl(k) = num2(maxpos(1));
        
        posiz{k} = maxpos;
        freq(k,:) = p;
        
        clear p y
    end
    
    f1 = 1:clust_number1;
    f2 = 1:clust_number2;
    
    for w=1:clust_number1
        maxx = max(max(freq(f1,f2)));
        indmax = find(freq(f1,f2)==maxx);
        indmax = indmax(1);
        clL2 = ceil(indmax/length(f1));
        clL1 = indmax-length(f1)*floor(indmax/length(f1));
        if clL1==0
            clL1 = length(f1);
        end
        indL2 = find(L2==num2(f2(clL2)));
        new_vector(indL2) = num1(f1(clL1));
        nuovocl(w) = num1(f1(clL1));
        usedcl(w) = num2(f2(clL2));
        f1(clL1) = [];
        f2(clL2) = [];
        
        clear indL2
    end
    
    missinglabel = num2(not(ismember(num2,nuovocl)));
    missingind = num2(not(ismember(num2,usedcl)));
    intersec = missinglabel(ismember(missinglabel,missingind));
    if not(isempty(intersec))
        missinglabel(ismember(missinglabel,intersec)) = [];
        missingind(ismember(missingind,intersec)) = [];
        for jj=1:length(intersec)
            indL2 = find(L2==intersec(jj));
            new_vector(indL2) = intersec(jj);
            clear indL2
        end
    end
    
    for j=1:length(missinglabel)
        if not(isempty(missingind)) && length(find(new_vector==0))>0
            indL2 = find(L2==missingind(j));
            new_vector(indL2) = missinglabel(j);
            clear indL2
        end
    end
    
end




end

