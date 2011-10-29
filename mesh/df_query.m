function [id] = df_query( DF, p)
    p = p-DF.MIN;
    p = p./DF.S;
    p = round(p)+1;
    
    id = DF.N(p(1),p(2),p(3));
end

