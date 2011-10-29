function [DF] = df_build(V,F,C)
    
    if ~exist('C','var')
        C=[10;10;10];
    end
    
    MIN = min(V,[],1);
    MAX = max(V,[],1);

    S = (MAX-MIN)./C';
    
    D = zeros(C'+1);
    N = zeros(C'+1);
    
    for x=1:size(D,1)
        progressbar(x,size(D,1));
        for y=1:size(D,2)
            for z=1:size(D,3)
                i = [x,y,z];
                p = MIN + (i-1).*S;
                [D(x,y,z) N(x,y,z)] = min(normrow(repmat(p,size(V,1),1)-V));
            end
        end
    end
    
    
    % Prepare DF struct
    DF.MIN = MIN;
    DF.MAX = MAX;
    DF.D = D;
    DF.N = N;
    DF.C = C;
    DF.S = S;
end