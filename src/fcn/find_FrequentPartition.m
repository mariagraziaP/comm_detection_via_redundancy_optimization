function [freqPart, pos] = find_FrequentPartition(partitions)
% partition = [Nodes x (NLayers) x NIterations]
% freqPart = [Nodes x (NLayers)]

if length(size(partitions))<3 % single layer community detection
    curr_ris = partitions;
    for i=1:size(partitions, 2)
%         sprintf('it = %g', i)
        sim_vi = zeros(size(partitions,2), 1);
        for j=1:size(partitions, 2)
            sim_vi(j) = partition_distance(...
                    curr_ris(:,i),curr_ris(:,j));
        end
        sumsim_vi(i) = sum(sim_vi);
    end
    
    freqPart = curr_ris(:,find(sumsim_vi==min(sumsim_vi),1));
    pos = find(sumsim_vi==min(sumsim_vi),1);

else % multilayer community detection
    for l=1:size(partitions, 2)
        curr_ris = squeeze(partitions(:,l,:));
        for i=1:size(partitions, 3)
%             sprintf('l = %g, it = %g', l, i)
            sim_vi = zeros(size(partitions,3), 1);
            for j=1:size(partitions, 3)
                sim_vi(j) = partition_distance(...
                    curr_ris(:,i),curr_ris(:,j));
            end
            sumsim_vi(i,l) = sum(sim_vi);
        end
        freqPart(:,l) = curr_ris(:,find(sumsim_vi(:,l)==min(sumsim_vi(:,l)),1));
    end
    
end
end
