% batch test    
clear;
name_list = {'constant','mountain','road2'};

for i=1:length(name_list)
    name = name_list{i};
    % test_landmark_error;
    test_distance_approximation;
end