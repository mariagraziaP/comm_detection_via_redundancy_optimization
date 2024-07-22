clearvars
close all

%% set path

addpath(genpath('...\fcn'))
addpath(genpath('...\ext'))

%% set directories 

data_dir = '...\mat';
out_dir = '...\out';
comm_dir = '...\out\optimized_comm';

%% load data and parameters

load(fullfile(data_dir, 'grandaverage_HCP'))
load(fullfile(data_dir, 'yeo7_200'));

randorder = randperm(N);
FCr = FC(randorder,randorder);  % randomize node order
yeo7r = yeo7(randorder);

N = size(FCr,1);

%% compute TSE curve

si = 10000;
[C, Ilvl, Ilvl_max] = calcCalt_Nint(FC, si);
[O, TC, D] = calcO_logdet(FC);

save(fullfile(out_dir, 'TSE_HCP_curve'), 'C', 'Ilvl', 'Ilvl_max', 'O', 'TC', 'D')


%% set optimization parameters

iter = 100;         % number of independent attempts (should be multiple of number of available PC cores)
riter = 1000000;    % number of iterations per attempt

%% annealing variables (vary hfrac as needed btw 10 and 100)

H = riter; hfrac = 10; Texp = 1-hfrac/H; T0 = 1;
flag = 0;  % display cost
Tall = T0.*Texp.^[1:1:H];

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% MAIN OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% loop over numclust

ncall = [];
tic

for numclust = 2:12
    
    ncall = [ncall numclust];
    disp(num2str(numclust))
        
    % optimization of TC within modules

    cost_out = zeros(1,iter);
    C_out = zeros(N,iter);
    
    parfor it= 1:iter
        
        % plant a seed and compute inital cost
        C0 = randi([1 numclust],N,1);
        
        mod_size = histcounts(C0, [0.5:1:numclust+0.5]);        
        TCmod = zeros(1,numclust);
        for m=1:numclust
            ff = find(C0==m);
            %[~, TCmod(m), ~] = calcO_logdet(FC(ff,ff));
            TCmod(m) = calcI_alt(FCr(ff,ff));
        end        
        cost0 = -mean(TCmod - Ilvl(mod_size));
        
        % initialize optimization
        h = 0; hcnt = 0;

        Copt = C0;
        costopt = cost0;
        
        inicost = cost0;
        maxcost = inicost;
        hicost = maxcost;
        
        Cmax = C0;
        
        while h<H
        
            h = h+1; hcnt = hcnt+1;
        
            %tt(h) = maxcost;
            
            % current temperature
            Tc = T0*Texp^h;
        
            if ((mod(h,H/10)==0)&&(flag==1))
                disp(['at step ',num2str(h),' - elapsed time ',num2str(toc),' - lowest cost = ',num2str(mincost)]);
            end
        
            rp = randperm(N);
            %node = rp(1:ceil(abs(randn)));
            node = rp(1);       % mutate one node
            Ctest = Cmax;
            
            for k=1:length(node)
                pickfrom = setdiff(1:numclust,Ctest(node(k)));
                rr = randi(numclust-1,1);
                Ctest(node(k)) = pickfrom(rr);
            end
        
            mod_size = histcounts(Ctest, [0.5:1:numclust+0.5]);
            TCmod = zeros(1,numclust);
            for m=1:numclust
                ff = find(Ctest==m);
                %[~, TCmod(m), ~] = calcO_logdet(FC(ff,ff));
                TCmod(m) = calcI_alt(FCr(ff,ff));
            end
            costtest = -mean(TCmod - Ilvl(mod_size)); 
            costnew = costtest;
        
            % annealing
            randcon = rand < exp(-(costnew-hicost)/Tc);
            %disp([num2str(costnew),' ',num2str(lowcost),' ',num2str(mincost),' ',num2str(costnew<lowcost),' ',num2str(randcon)]);
            if (costnew < hicost) || randcon
                Cmax = Ctest;
                hicost = costnew;
                % is this a new absolute best?
                if (hicost < maxcost)
                    Cmaxglobal = Cmax;
                    maxcost = hicost;
                    hcnt = 0;
                end
            end
        end
        
        C_out(:,it) = Cmaxglobal;
        cost_out(it) = -maxcost;
        
        disp([num2str(numclust),' | ',num2str(-maxcost)]);
                
    % end parfor
    end 
    
    % peak cost
    [max1,max2] = max(cost_out);

    % how similar are the 'iter' community vectors to each other?
    if iter>1
        [v, m] = partition_distance(C_out);
        utri = find(triu(ones(iter),1));
        comm_sim(numclust) = mean(m(utri));
    end
    
    % how similar are the 'iter' community vectors to 'yeo7'?
    [v2, m2] = partition_distance(C_out,yeo7r);
    sim2yeo(numclust) = mean(m2);
    sim2yeo_max(numclust) = m2(max2);
    
    disp(['num mods = ',num2str(numclust),'  | mean cost = ',num2str(mean(cost_out)),'  | max cost = ',num2str(max(cost_out))]);


    %% compute o-info between communities
    
    samp = 100;
    
    for i=1:iter
        
        CI = C_out(:,i);
        
        for s=1:samp
            for nc=1:numclust
                ff = find(CI==nc);
                rp = randperm(length(ff));
                nod(nc) = ff(rp(1));
            end
            [O(s), I(s), D(s)] = calcO_logdet(FCr(nod,nod));
            clear nod
        end
        
        Obc(i) = mean(O); Ibc(i) = mean(I); Dbc(i) = mean(D);
        
    end
    cc(numclust) = corr(cost_out',Obc','type','spearman');
    
    % collect all output
    C_out_all(:,:,numclust) = C_out;
    cost_out_all(:,numclust) = cost_out;
    Obc_all(:,numclust) = Obc;
    Ibc_all(:,numclust) = Ibc;
    Dbc_all(:,numclust) = Dbc;    
    Ibc2_all(:,numclust) = Ibc-Ilvl(numclust);
    
    % ratio TC over mean(O) at peak TC
    TCO_rat(numclust) = max1/Obc_all(max2,numclust);
    TCI_rat(numclust) = max1/Ibc_all(max2,numclust);
    
    % co-classification matrix
    Aall = agreement(C_out)./iter;
    % analytic null
    anull = 0;
    for cnum = 1:size(C_out,2)
        anull = anull+sum((histcounts(C_out(:,cnum))./N).*((histcounts(C_out(:,cnum))-1)./(N-1)));
    end
    A = Aall - anull/iter;
    % consensus clustering (specific for each value of nclust)
    tau = 0;
    ciu(:,numclust) = consensus_und(A,tau,10);   
    
end

tempo = toc;


%% save

save(fullfile(comm_dir, sprintf('anneal_comm_%giter', iter)), 'randorder', 'C_out_all',...
    'cost_out_all', 'ciu', 'Obc_all', 'Ibc_all', 'Dbc_all', 'Ibc2_all',...
    'comm_sim', 'sim2yeo', 'sim2yeo_max', 'TCO_rat', 'TCI_rat', 'tempo',...
    'H', 'T0', 'hfrac', 'Texp', 'riter', 'iter', 'samp')


