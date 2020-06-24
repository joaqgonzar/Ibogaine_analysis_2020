%% Clsuter-based multiple comparisons correction
% 
% This program implements a multiple comparisons correction by means of
% finding significant frequency clusters (by means of a paired t-test)
% and then applying a permutation test to find a statistical threshold.
% Thus, instead of pre defining the requency bands of an EEG, these emerge
% directly from data as the clusters of frequencies that change between the  
% 2 conditions which are above chance.  
% 
% The implementation is made as following: an initial paired t-test is
% carried on a frequency by frequency manner, to obtain the significant
% frequencies. To control the multiple comparisons problem, all the
% significant frequencies adjacent to one another are grouped into
% clusters. A new statistic per cluster is obtained by summing the 
% individual t-stats from each frequency inside a cluster. A null
% distribution for the clusters statistic is constructed by randomizing 
% the labels of each animal (participant). The p_value for each cluster is
% then find comparing the cluster statistic against the null distribution.
% 
% 
%
% The inputs are the following: 
% 
% power1: the power spectrums in the first condition condition (power1),
% which should be formatted as power1(animals, frequencies). 
% power2: the power spectrums in the second condition condition (power2),
% which should be formatted as power1(animals, frequencies).
% tail: the tail of the test to be employed, either 'both', 'right' or 'left'
% minimum cluster size: minimum number of continous frequencies to be taken
% as a cluster 
% num_permutations: number of randomizations to be made in the correction
%
% The outputs are the following: 
% 
% freq_clusters: frequency clusters detected
% p_value_clusters: the p_values of each freq_cluster
% stat_cluster: the statistic of each cluster (sum of the t-stats whithin 
% the cluster)
%
% References: Gonzalez, et al, 2020c bioRxiv.  Maris & Oostenveld, 2007,
% Journal of Neuroscience Methods
% 
% Joaquin Gonzalez, 2020, Laboratorio de Neurobiologia del Sueno, Facultad 
% de Medicina, Universidad de la Republica. email: joaqgonzar@gmail.com



function [freq_clusters,p_value_clusters,stat_cluster] = Cluster_Permutation_Correction_xfreq(power1,power2,tail,min_clast_size,num_permutations)
 
    
    % Start by computing the standard t-test in each frequency
    [h,p,ci,stat] = ttest(power1,power2,'tail',tail); % ttest con la cola elegida
    h = squeeze(h); 
    p = squeeze(p); 
    estadistico = squeeze(stat.tstat);
    ranking = sort(p);

    % CLUSTER BASED PERMUTATION CORRECTION
    significant = find(h>0);
    differences = diff(significant);
    clusters = find(differences>1);
    
    % single cluster case
    if sum(significant)>0 & sum(differences>1)==0
        clusters = 1;
    else 
    end
    
    % multiple clusters extraction
    if length(clusters) == 0;
        display('No clusters found')
        freq_clusters = [];
        p_value_clusters = [];
        stat_cluster = [];
        return
    else 
    end
    
    index = 1;
    for i = 1:length(clusters)
            clus{i} = num2cell(significant(index:clusters(i)));
            index = clusters(i)+1;
        if i == length(clusters);
            clus{i+1} = num2cell(significant(index:end));
        else
        end
    end

    %
    % control for cluster size
    index_clusters = 1;
    for i = 1:length(clus)
        if size(clus{i},2) > min_clast_size;
            freq_clusters{index_clusters} = clus{i};
            index_clusters = index_clusters+1;
        else 
        end
    end
    
    if index_clusters > 1
    
    number_cluster = freq_clusters(:);
    number_cluster = sum(~cellfun(@isempty,number_cluster));
    
    for i = 1:number_cluster;
            clus = freq_clusters{i};
            tmat = cell2mat(clus);
            stat_cluster(i) = num2cell(sum(estadistico(tmat)));
    end

    % VERSION ALEATORIZADA, clear 
    
    for j = 1:num_permutations;
        data = cat(1,power1,power2);

        data_rand = randsample([1:size(data,1)],size(data,1),true);
        data1 = data(data_rand(1:size(data,1)/2),:);
        data2 = data(data_rand(size(data,1)/2 + 1:end),:);
        %
        [h,p,ci,stat] = ttest(data1,data2,'tail','both'); % ttest con la cola elegida
        h = squeeze(h);
        p = squeeze(p); 
        estadistico_t = squeeze(stat.tstat);

        % CLUSTER BASED PERMUTATION CORRECTION
        
        significant = find(h>0);
        differences = diff(significant);
        clusters = find(differences>1);

        
        index = 1;
        for i = 1:length(clusters)
                clus{i} = num2cell(significant(index:clusters(i)));
                index = clusters(i)+1;
            if i == length(clusters);
                clus{i+1} = num2cell(significant(index:end));
            else
            end
        end

        % control por tama??o de cluster
        index_clusters = 1;
        for i = 1:length(clus)
            if size(clus{i},2) >= min_clast_size;
                clusters_final{index_clusters} = clus{i};
                index_clusters = index_clusters+1;
            else 
            end
        end
        
        if index_clusters > 1
        for i = 1:length(clusters_final);
                clus = clusters_final{i};
                tmat = cell2mat(clus);
                stat_cluster_random(i,j) = sum(estadistico(tmat));
        end
        else 
        end
    end

    % histograma de la distribucion nula
    stat_cluster_random = reshape(stat_cluster_random,[1 size(stat_cluster_random,1)*size(stat_cluster_random,2)]);
    %
    for i = 1:number_cluster
        pr = sum(stat_cluster_random >= cell2mat(stat_cluster(i)))/length(stat_cluster_random);
        pi = sum(stat_cluster_random <= cell2mat(stat_cluster(i)))/length(stat_cluster_random);
        
        if strcmp(tail,'right')
            p_value_clusters(i) = num2cell(pr); 
        else 
            p_value_clusters(i) = num2cell(2*min(pr,pi));
        end 
           
        if cell2mat(p_value_clusters(i)) == 0
            p_value_clusters(i) = num2cell(1/length(stat_cluster));
        else
        end
    end
    else
        display('No clusters were found')
        freq_clusters = [];
        p_value_clusters = [];
        stat_cluster = [];
        
    end
    %p_valor_clusters
    

end