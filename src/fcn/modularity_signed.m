function Q = modularity_signed(adj, com, gamma)

% adj = adjacency matrix, undirected signed
% com = community membership

N = length(adj);

Apos = adj.*(adj>0);
Aneg = -adj.*(adj<0);

% positive weights

Kp = sum(Apos);                               %degree
mp = sum(Kp);                               %number of edges (each undirected edge is counted twice)
% Bp = Apos - gamma*(Kp.'*Kp)/mp;                    %modularity matrix
Bp = Apos - gamma;                    %modularity matrix, uniform null model
sp = com(:,ones(1,N));                      %compute modularity
Qp =~ (sp-sp.').*Bp/mp;
Qp = sum(Qp(:));

% negative weights

Kn = sum(Aneg);                               
mn = sum(Kn)+sum(Kp);                              
% Bn = Aneg - gamma*(Kn.'*Kn)/mn;        
Bn = Aneg - gamma;                  
sn = com(:,ones(1,N));                      
Qn =~ (sn-sn.').*Bn/mn;
Qn = sum(Qn(:));

% sum

Q = Qp - Qn;

end





