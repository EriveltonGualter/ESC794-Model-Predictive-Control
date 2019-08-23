% function animate(t,x,P)
function animate(car1, car2, car3, car4, t)

    targetFrameRate = 10;

    % Animation call-back variables:
    IS_PAUSED = false;
    SPEED = 1;
    QUIT = false;
    START_TIME = t(1);
    SIM_TIME = START_TIME;
    SIMULATION_COMPLETE = false;

    fig = figure(2);
    % UIControls
    slr = uicontrol(fig,'Style','slider');
    slr.Position = [20 75 70 20];
    slr.Callback = @selection;
    set(slr, 'min', 0.05);
    set(slr, 'max', 2);
    set(slr, 'Value', 1); 
    set(slr, 'SliderStep', [0.05, 0.1]);

    btn1 = uicontrol(fig,'Style','pushbutton'); 
    btn1.Position = [20 100 70 20];
    btn1.String = {'Restart'};
    btn1.Callback = @(src, evetn)restart_sim(src, evetn, IS_PAUSED, car1, car2, car3, car4, t);

    btn2 = uicontrol(fig,'Style','pushbutton'); 
    btn2.Position = [20 125 70 20];
    btn2.String = {'Pause'};
    btn2.Callback = @(src, evetn)pause_sim(src, evetn, IS_PAUSED);

    btn3 = uicontrol(fig,'Style','pushbutton'); 
    btn3.Position = [20 150 70 20];
    btn3.String = {'Play'};
    btn3.Callback = @(src, evetn)play_sim(src, evetn, IS_PAUSED, car1, car2, car3, car4, t);

    btn4 = uicontrol(fig,'Style','pushbutton'); 
    btn4.Position = [20 175 70 20];
    btn4.String = {'Quit'};
    btn4.Callback = @quit_sim;

    txt = uicontrol(fig,'Style','text'); 
    txt.Position = [20 50 70 20];
    txt.String = {'1x'};

    sim_an(car1, car2, car3, car4, t);

    function quit_sim(src,evetn,ip)
        QUIT = true;
    end

    function play_sim(src,evetn,ip, car1, car2, car3, car4, t)
        IS_PAUSED = false;
        if SIMULATION_COMPLETE
            SIM_TIME = START_TIME;
            sim_an(car1, car2, car3, car4, t);
        end
    end

    function pause_sim(src,evetn,ip)
        IS_PAUSED = true;
    end

    function restart_sim(src,event, ip, car1, car2, car3, car4, t)
        SIM_TIME = START_TIME;
        if SIMULATION_COMPLETE
            sim_an(car1, car2, car3, car4, t);
        end
    end

    function selection(src,event)
        SPEED = slr.Value;
        str = horzcat('x',num2str(SPEED));
        txt.String = str;
    end

    function sim_an(car1, car2, car3, car4, t)
        tic;    %Start a timer
        timeBuffer(1:3) = toc;
        while SIM_TIME < t(end);
            SIMULATION_COMPLETE = false;
            pcar1 = interp1(t',car1',SIM_TIME,'linear','extrap')';
            pcar2 = interp1(t',car2',SIM_TIME,'linear','extrap')';
            pcar3 = interp1(t',car3',SIM_TIME,'linear','extrap')';
            pcar4 = interp1(t',car4',SIM_TIME,'linear','extrap')';

            cla
            drawStreet(SIM_TIME)   
            drawCar(pcar1(1), pcar1(2),     [1 0 0]); 
            drawCar(pcar2(1), pcar2(2),     [0 0 1]); 
            drawCar(pcar3(1), pcar3(2),     [1 1 0]); 
            drawCar(pcar4(1), pcar4(2),     [1 0 1]); 

            drawnow;

            %Set up targets for timing
            dtReal = 0.5*(timeBuffer(1) - timeBuffer(3));
            if IS_PAUSED
                dtSim = 0;
            else
                dtSim = SPEED*dtReal;
            end
            SIM_TIME = SIM_TIME + dtSim;

            %Record the frame rate:
            timeBuffer(3) = timeBuffer(2);
            timeBuffer(2) = timeBuffer(1);
            timeBuffer(1) = toc;

            % Check exit conditions:
            if QUIT
                close all
                break
            end

        end
        SIMULATION_COMPLETE = true;
    end

end %animate.m

function drawCar(x,y,c)    
    if abs(x) == 1
        a = 1; b = 1.2;
    else
        b = 1; a = 1.2;
    end

    rectangle(...
            'Position',[x-0.5*a, y-0.5*b, a, b],...
            'Curvature',0.1,...
            'LineWidth',2,...
            'FaceColor',0.8*c,...
            'EdgeColor',c);  %Darker version of color
        
    axis equal; axis([-12 12 -12 12]); set(gca,'xtick',[]); set(gca,'ytick',[]);
end

function drawStreet(time)
    title(sprintf('t = %2.2f%  s',time));
    rectangle('Position',[-10 2 8 8],   'Curvature', 0.1, 'LineWidth',2, 'FaceColor',0.8*[0 1 0], 'EdgeColor',0.5*[0 1 0]); 
    rectangle('Position',[-10 -10 8 8], 'Curvature', 0.1, 'LineWidth',2, 'FaceColor',0.8*[0 1 0], 'EdgeColor',0.5*[0 1 0]);  
    rectangle('Position',[2 -10 8 8],   'Curvature', 0.1, 'LineWidth',2, 'FaceColor',0.8*[0 1 0], 'EdgeColor',0.5*[0 1 0]); 
    rectangle('Position',[2 2 8 8],     'Curvature', 0.1, 'LineWidth',2, 'FaceColor',0.8*[0 1 0], 'EdgeColor',0.5*[0 1 0]); 
    line([-9 -3],[0 0],'Color','white','LineStyle','--','LineWidth',2)
    line([ 9  3],[0 0],'Color','white','LineStyle','--','LineWidth',2)
    line([0 0], [-9 -3],'Color','white','LineStyle','--','LineWidth',2)
    line([0 0], [ 9  3],'Color','white','LineStyle','--','LineWidth',2)
    axis([-12 12 -12 12])
    axis equal
end


