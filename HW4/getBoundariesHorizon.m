function Xbd = getBoundariesHorizon(Ad, Bd, step, N)
% %UNTITLED2 Summary of this function goes here
% %   Detailed explanation goes here
  
    % Initialize
    uarray = linspace(-1,1,step);
    ii = 1;
    
    for xi=linspace(-1,1,step)
        for yi=linspace(-1,1,step)
            X0(ii,:) = [xi yi];
            ii = ii+1;
        end
    end
    
    % Create all possible permutations (with repetition) for u control 
    u = linspace(-1,1,8);
    C = cell(8, 1);             % Preallocate a cell array
    [C{:}] = ndgrid(u);         % Create K grids of values
    uall = cellfun(@(x){x(:)}, C); % Convert grids to column vectors
    uall = [uall{:}];   

    area = [];
    for i=1:length(X0)
        X = X0(i,:)';
        
        j = 1;
        flag = 1;
        Xp = [];
        while flag  % Flag: Reachs the length of uarray
            Xc = X;                         % Current X
            for k=1:N
                Xv = Ad*Xc + Bd*uall(j,k);   % Next X
                Xc = Xv;
                Xp = [Xp Xv];
                if (abs(Xv(1))>1 || abs(Xv(2))>1)
                    break
                end
                if (abs(Xv(1))<=1 && abs(Xv(2))<=1) && norm(Xv)<=0.5
                    area = [area X];
                    flag = 0;
                    break
                end

            end
            if j == length(uarray)
                flag = 0;
            else
                j = j+1;
            end
        end        
    end
    
    bd = boundary(area(1,:)', area(2,:)');
    Xbd = [area(1,bd) area(2,bd)];
    Xbd = [area(1,bd); area(2,bd)];

    % Plot
    figure; hold on;
    rectangle('Position',[-1 -1 2 2],'EdgeColor','r', 'LineWidth',2)
    rectangle('Position',[-0.5 -0.5 1 1],'EdgeColor','b', 'LineWidth',2,'Curvature',[1 1])
    axis equal
    s = scatter(area(1,:),area(2,:),'filled'); s.MarkerFaceColor = 'r'; hold off
    alpha(s,.01); hold on;
    plot(area(1,bd),area(2,bd), 'k','LineWidth', 3);
    axis([-1.5 3 -1.5 1.5]);  axis equal
    text(1.2,0.25,'A = ');    text(1.4,0, num2str(round(Ad,2)))
    text(1.2,-.35,'B = ');    text(1.4,-0.6, num2str(round(Bd,2)))
    Nstr = horzcat('N = ', num2str(round(N,2))); text(1.2,-0.9, Nstr)
    
end


