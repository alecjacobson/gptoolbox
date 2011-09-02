function vectarrow(P,V)
%   P: list of points #x3
%   V: list of vectors #x3
%   Rentian Xiong 4-18-05

  for index = 1:size(P,1)
    hold on;
    p0 = P(index,:);
    p1 = P(index,:) + V(index,:);

    if max(size(p0))==3
        if max(size(p1))==3
            x0 = p0(1);
            y0 = p0(2);
            z0 = p0(3);
            x1 = p1(1);
            y1 = p1(2);
            z1 = p1(3);
            plot3([x0;x1],[y0;y1],[z0;z1]);   % Draw a line between p0 and p1
            
            p = p1-p0;
            alpha = 0.1;  % Size of arrow head relative to the length of the vector
            beta = 0.1;  % Width of the base of the arrow head relative to the length
            
            hu = [x1-alpha*(p(1)+beta*(p(2)+eps)); x1; x1-alpha*(p(1)-beta*(p(2)+eps))];
            hv = [y1-alpha*(p(2)-beta*(p(1)+eps)); y1; y1-alpha*(p(2)+beta*(p(1)+eps))];
            hw = [z1-alpha*p(3);z1;z1-alpha*p(3)];
            
            hold on
            plot3(hu(:),hv(:),hw(:))  % Plot arrow head
            grid on
            xlabel('x')
            ylabel('y')
            zlabel('z')
            hold off
        else
            error('p0 and p1 must have the same dimension')
        end
    elseif max(size(p0))==2
        if max(size(p1))==2
            x0 = p0(1);
            y0 = p0(2);
            x1 = p1(1);
            y1 = p1(2);
            plot([x0;x1],[y0;y1]);   % Draw a line between p0 and p1
            
            p = p1-p0;
            alpha = 0.1;  % Size of arrow head relative to the length of the vector
            beta = 0.1;  % Width of the base of the arrow head relative to the length
            
            hu = [x1-alpha*(p(1)+beta*(p(2)+eps)); x1; x1-alpha*(p(1)-beta*(p(2)+eps))];
            hv = [y1-alpha*(p(2)-beta*(p(1)+eps)); y1; y1-alpha*(p(2)+beta*(p(1)+eps))];
            
            hold on
            plot(hu(:),hv(:))  % Plot arrow head
            grid on
            xlabel('x')
            ylabel('y')
            hold off
        else
            error('p0 and p1 must have the same dimension')
        end
    else
        error('this function only accepts 2D or 3D vector')
    end
    hold off;
  end
