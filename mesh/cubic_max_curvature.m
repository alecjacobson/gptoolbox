function max_k = cubic_max_curvature(C,t0,t1)
% Max |curvature| of cubic Bézier C on [t0,t1]
% C is 4x2

assert(size(C,1)==4 && size(C,2)==2);

%% ------------------------------------------------------------
% Step 1: find inflection points in (t0,t1)
%% ------------------------------------------------------------

P = C;
Cidx = 1:4;

[I,Tinf] = spline_inflection_points(P,Cidx);

% keep strictly inside interval
Tinf = Tinf(Tinf>t0 & Tinf<t1);
Tinf = sort(Tinf);

%% ------------------------------------------------------------
% Step 2: build parameter breakpoints
%% ------------------------------------------------------------

Tbreak = [t0; Tinf; t1];
max_k = 0;

%% ------------------------------------------------------------
% Step 3: process each sub-interval
%% ------------------------------------------------------------

for i = 1:length(Tbreak)-1

    ta = Tbreak(i);
    tb = Tbreak(i+1);

    % Map [ta,tb] to [0,1]
    if ta == 0 && tb == 1
        Csub = C;
    else
        % First split at tb
        [Cleft,~] = cubic_split(C,tb);

        % Then split left piece at ta/tb
        if ta > 0
            s = ta / tb;
            [~,Csub] = cubic_split(Cleft,s);
        else
            Csub = Cleft;
        end
    end

    % Now compute max curvature on full [0,1]
    k = cubic_max_curvature_single_interval(Csub);

    max_k = max(max_k,k);
end

end
