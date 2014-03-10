name_list = {'bunny','elephant-50kv','david50kf','david_head','hand'};

for nstart=[1 10 20 50 100]
    for iname = 1:length(name_list)
        name = name_list{iname};
        disp(['----> Processing mesh ' name ' - ' num2str(nstart) '.']);
        test_propagation_mesh;
    end
end