%% Clsuter-based multiple comparisons correction
% 
% This program provides an example of usage to the funcion 
% Cluster_Permutation_Correction_xfreq by generating a random power 
% spectrums and then introdcing a power increase in 2 different set of 
% frequencies.
% 
% Joaquin Gonzalez, 2020, Laboratorio de Neurobiologia del Sueno, Facultad 
% de Medicina, Universidad de la Republica. email: joaqgonzar@gmail.com

power1 = randn(6,1000); % random data
power2 = randn(6,1000); % random data

power2(:,100:200) = 10; % introduce a power increase between 100 and 200 'Hz' 
power2(:,500:550) = 10; % introduce a power increase between 500 and 550 'Hz' 

tail = 'both';
min_clast_size = 10;
num_permutations = 1000;

[freq_clusters,p_value_clusters,stat_cluster] = Cluster_Permutation_Correction_xfreq(power1,power2,tail,min_clast_size,num_permutations);

freq_clusters{:}
p_value_clusters