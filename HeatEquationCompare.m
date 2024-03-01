function [] = HeatEquationCompare()
    % HeatEquationCompare(m,n,TT,dt,dx,D)
    % Initial variables, the space are splitted into mxn grid, (m+1)x(n+1) points.
    % Note: dx = dy, so we can use dx instead of dy.
    %{
    m = 20;
    n = 20;
    TT = 1;
    dt = 0.01;
    dx = 0.1;
    dy = 0.1;
    D = 0.1;
    %}

    % Define the number of files
    numFiles = 10;

    % Loop over the file names
    for i = 1:numFiles
        % Create the file name
        fileName = sprintf('testcases/%d.mat', i);
        
        % Load the .mat file
        load(fileName);

        errMax = 0.5;
        timeStep = TT/dt;

        %{
        w = 0; % Parameter of SOR solver
        T_mid_max = 0;

        % Find w with maximum center point for SOR solver 
        for w_range=0.01:0.01:1.99
            T_mid_tmp = 0;
            % Initialize T
            T = zeros(m+1,n+1);
            FD = zeros(m+1,n+1);
            T(1,2:end-1) = 900;   % Top
            T(end,2:end-1) = 600; % Bottom
            T(2:end-1,1) = 400;   % Left
            T(2:end-1,end) = 800; % Right
            T(1,1) = 650; T(1,end) = 850; T(end,end) = 700; T(end,1) = 500; % Corner nodes
            % First time step
            % Loop until the error of T is small enough (< errMax)
            err = 1;
            while(err >= errMax)
                T_old = T;
                for x = 2:m
                    for y = 2:n
                        T(x,y) = (1-w)*T_old(x,y) + w*1/4*(T(x-1,y)+T_old(x+1,y)+T(x,y-1)+T_old(x,y+1)); % SOR formula
                    end
                end
                err = max(abs(T-T_old),[],'all');
            end
        
            % Calculate FD
            for x = 2:m
                for y = 2:n
                    FD(x,y) = (T(x-1,y)+T(x+1,y)+T(x,y-1)+T(x,y+1)-4*T(x,y))/(dx*dx);
                end
            end
    
            % For other time steps
            for t = 2:timeStep+1
                % Calculate new T
                for x = 2:m
                    for y = 2:n
                        T(x,y) = FD(x,y)*dt*D + T(x,y); % Forward Euler
                    end
                end
        
                % Loop until the error of T is small enough (< errMax)
                err = 1;
                while(err >= errMax)
                    T_old = T;
                    for x = 2:m
                        for y = 2:n
                            T(x,y) = (1-w)*T_old(x,y) + w*1/4*(T(x-1,y)+T_old(x+1,y)+T(x,y-1)+T_old(x,y+1)); % SOR formula
                        end
                    end
                    err = max(abs(T-T_old),[],'all');
                end
            
                % Calculate FD
                for x = 2:m
                    for y = 2:n
                        FD(x,y) = (T(x-1,y)+T(x+1,y)+T(x,y-1)+T(x,y+1)-4*T(x,y))/(dx*dx);
                    end
                end
            end
            T_mid_tmp = T(m/2+1,n/2+1);
            if T_mid_tmp > T_mid_max
                w = w_range;
                T_mid_max = T_mid_tmp;
            end
        end
        %}
        w = 1.26; % Parameter of SOR solver

    
        % Initialize T
        T = zeros(m+1,n+1);
        T(1,2:end-1) = 900;   % Top
        T(end,2:end-1) = 600; % Bottom
        T(2:end-1,1) = 400;   % Left
        T(2:end-1,end) = 800; % Right
        T(1,1) = 650; T(1,end) = 850; T(end,end) = 700; T(end,1) = 500; % Corner nodes
    
        % Initialize FD
        FD = zeros(m+1,n+1);
        T_mid = zeros(timeStep+1,1);
    
        % Jacobi solver
        % First time step
        % Loop until the error of T is small enough (< errMax)
        cnt1 = zeros(timeStep+1,1);
        err = 1;
        while(err >= errMax)
            T_old = T;
            for x = 2:m
                for y = 2:n
                    T(x,y) = 1/4*(T_old(x-1,y)+T_old(x+1,y)+T_old(x,y-1)+T_old(x,y+1)); % Approximate formula for 2nd order derivative
                end
            end
            err = max(abs(T-T_old),[],'all');
            cnt1(1)=cnt1(1)+1;
        end
    
        % Calculate FD
        for x = 2:m
            for y = 2:n
                FD(x,y) = (T(x-1,y)+T(x+1,y)+T(x,y-1)+T(x,y+1)-4*T(x,y))/(dx*dx);
            end
        end
        T_mid(1) = T(m/2+1,n/2+1);
    
        % Plot filled contours with color and contour lines with labels
        contourf(T);
        colorbar;
        colormap(jet);
        
        % Add contour lines with labels
        hold on;
        contour(T, 'LineColor', 'k', 'ShowText', 'on');
        
        % Customize the plot
        xlabel('x');
        ylabel('y');
        title('Jacobi Solver');
    
        % Capture the current frame
        frame = getframe(gcf);
    
        % Convert the frame to an indexed image
        im = frame2im(frame);
        
        % Convert the indexed image to a colormap
        [imind, cm] = rgb2ind(im, 256);

        fileName = sprintf('outputs/%d/JacobiSolver.gif', i);
    
        imwrite(imind, cm, fileName, 'gif', 'Loopcount', inf);
        
        % Pause between frames (optional)
        % pause(0.01);
    
        clf;
    
        % For other time steps
        for t = 2:timeStep+1
            cnt1(t)=cnt1(t-1);
            % Calculate new T
            for x = 2:m
                for y = 2:n
                    T(x,y) = FD(x,y)*dt*D + T(x,y); % Forward Euler
                end
            end
    
            % Loop until the error of T is small enough (< errMax)
            err = 1;
            while(err >= errMax)
                T_old = T;
                for x = 2:m
                    for y = 2:n
                        T(x,y) = 1/4*(T_old(x-1,y)+T_old(x+1,y)+T_old(x,y-1)+T_old(x,y+1)); % Approximate formula for 2nd order derivative
                    end
                end
                err = max(abs(T-T_old),[],'all');
                cnt1(t)=cnt1(t)+1;
            end
        
            % Calculate FD
            for x = 2:m
                for y = 2:n
                    FD(x,y) = (T(x-1,y)+T(x+1,y)+T(x,y-1)+T(x,y+1)-4*T(x,y))/(dx*dx);
                end
            end
            T_mid(t) = T(m/2+1,n/2+1);
    
            % Plot filled contours with color and contour lines with labels
            contourf(T);
            colorbar;
            colormap(jet);
            
            % Add contour lines with labels
            hold on;
            contour(T, 'LineColor', 'k', 'ShowText', 'on');
            
            % Customize the plot
            xlabel('x');
            ylabel('y');
            title('Jacobi Solver');
    
            % Capture the current frame
            frame = getframe(gcf);
    
            % Convert the frame to an indexed image
            im = frame2im(frame);
            
            % Convert the indexed image to a colormap
            [imind, cm] = rgb2ind(im, 256);

            fileName = sprintf('outputs/%d/JacobiSolver.gif', i);
    
            imwrite(imind, cm, fileName, 'gif', 'WriteMode', 'append');
        
            % Pause between frames (optional)
            % pause(0.01);
    
            % Clear the figure for the next frame
            clf;
        end
    
        T_mid_1 = T_mid;
    
        % Initialize T
        T = zeros(m+1,n+1);
        T(1,2:end-1) = 900;   % Top
        T(end,2:end-1) = 600; % Bottom
        T(2:end-1,1) = 400;   % Left
        T(2:end-1,end) = 800; % Right
        T(1,1) = 650; T(1,end) = 850; T(end,end) = 700; T(end,1) = 500; % Corner nodes
    
        % Gauss-Seidel solver
        % First time step
        % Loop until the error of T is small enough (< errMax)
        cnt2 = zeros(timeStep+1,1);
        err = 1;
        while(err >= errMax)
            T_old = T;
            for x = 2:m
                for y = 2:n
                    T(x,y) = 1/4*(T(x-1,y)+T_old(x+1,y)+T(x,y-1)+T_old(x,y+1)); % Gauss-Seidel solver formula
                end
            end
            err = max(abs(T-T_old),[],'all');
            cnt2(1) = cnt2(1)+1;
        end
    
        % Calculate FD
        for x = 2:m
            for y = 2:n
                FD(x,y) = (T(x-1,y)+T(x+1,y)+T(x,y-1)+T(x,y+1)-4*T(x,y))/(dx*dx);
            end
        end
        T_mid(1) = T(m/2+1,n/2+1);
    
        % Plot filled contours with color and contour lines with labels
        contourf(T);
        colorbar;
        colormap(jet);
        
        % Add contour lines with labels
        hold on;
        contour(T, 'LineColor', 'k', 'ShowText', 'on');
        
        % Customize the plot
        xlabel('x');
        ylabel('y');
        title('Gauss-Seidel Solver');
    
        % Capture the current frame
        frame = getframe(gcf);
    
        % Convert the frame to an indexed image
        im = frame2im(frame);
        
        % Convert the indexed image to a colormap
        [imind, cm] = rgb2ind(im, 256);
        
        fileName = sprintf('outputs/%d/GaussSeidelSolver.gif', i);
    
        imwrite(imind, cm, fileName, 'gif', 'Loopcount', inf);
    
        clf;
    
        % For other time steps
        for t = 2:timeStep+1
            cnt2(t) = cnt2(t-1);
            % Calculate new T
            for x = 2:m
                for y = 2:n
                    T(x,y) = FD(x,y)*dt*D + T(x,y); % Forward Euler
                end
            end
    
            % Loop until the error of T is small enough (< errMax)
            err = 1;
            while(err >= errMax)
                T_old = T;
                for x = 2:m
                    for y = 2:n
                        T(x,y) = 1/4*(T(x-1,y)+T_old(x+1,y)+T(x,y-1)+T_old(x,y+1)); % Gauss-Seidel solver formula
                    end
                end
                err = max(abs(T-T_old),[],'all');
                cnt2(t)=cnt2(t)+1;
            end
        
            % Calculate FD
            for x = 2:m
                for y = 2:n
                    FD(x,y) = (T(x-1,y)+T(x+1,y)+T(x,y-1)+T(x,y+1)-4*T(x,y))/(dx*dx);
                end
            end
            T_mid(t) = T(m/2+1,n/2+1);
    
            % Plot filled contours with color and contour lines with labels
            contourf(T);
            colorbar;
            colormap(jet);
            
            % Add contour lines with labels
            hold on;
            contour(T, 'LineColor', 'k', 'ShowText', 'on');
            
            % Customize the plot
            xlabel('x');
            ylabel('y');
            title('Gauss-Seidel Solver');
    
            % Capture the current frame
            frame = getframe(gcf);
    
            % Convert the frame to an indexed image
            im = frame2im(frame);
            
            % Convert the indexed image to a colormap
            [imind, cm] = rgb2ind(im, 256);

            fileName = sprintf('outputs/%d/GaussSeidelSolver.gif', i);
    
            imwrite(imind, cm, fileName, 'gif', 'WriteMode', 'append');
        
            % Pause between frames (optional)
            % pause(0.01);
    
            % Clear the figure for the next frame
            clf;
        end
    
        T_mid_2 = T_mid;
    
        % Initialize T
        T = zeros(m+1,n+1);
        T(1,2:end-1) = 900;   % Top
        T(end,2:end-1) = 600; % Bottom
        T(2:end-1,1) = 400;   % Left
        T(2:end-1,end) = 800; % Right
        T(1,1) = 650; T(1,end) = 850; T(end,end) = 700; T(end,1) = 500; % Corner nodes
    
        % SOR solver 
        % First time step
        % Loop until the error of T is small enough (< errMax)
        cnt3 = zeros(timeStep+1,1);
        err = 1;
        while(err >= errMax)
            T_old = T;
            for x = 2:m
                for y = 2:n
                    T(x,y) = (1-w)*T_old(x,y) + w*1/4*(T(x-1,y)+T_old(x+1,y)+T(x,y-1)+T_old(x,y+1)); % SOR formula
                end
            end
            err = max(abs(T-T_old),[],'all');
            cnt3(1) = cnt3(1)+1;
        end
    
        % Calculate FD
        for x = 2:m
            for y = 2:n
                FD(x,y) = (T(x-1,y)+T(x+1,y)+T(x,y-1)+T(x,y+1)-4*T(x,y))/(dx*dx);
            end
        end
        T_mid(1) = T(m/2+1,n/2+1);
    
        % Plot filled contours with color and contour lines with labels
        contourf(T);
        colorbar;
        colormap(jet);
        
        % Add contour lines with labels
        hold on;
        contour(T, 'LineColor', 'k', 'ShowText', 'on');
        
        % Customize the plot
        xlabel('x');
        ylabel('y');
        title('SOR Solver');
    
        % Capture the current frame
        frame = getframe(gcf);
    
        % Convert the frame to an indexed image
        im = frame2im(frame);
        
        % Convert the indexed image to a colormap
        [imind, cm] = rgb2ind(im, 256);

        fileName = sprintf('outputs/%d/SORSolver.gif', i);
    
        imwrite(imind, cm, fileName, 'gif', 'Loopcount', inf);
    
        clf;
    
        % For other time steps
        for t = 2:timeStep+1
            cnt3(t)=cnt3(t-1);
            % Calculate new T
            for x = 2:m
                for y = 2:n
                    T(x,y) = FD(x,y)*dt*D + T(x,y); % Forward Euler
                end
            end
    
            % Loop until the error of T is small enough (< errMax)
            err = 1;
            while(err >= errMax)
                T_old = T;
                for x = 2:m
                    for y = 2:n
                        T(x,y) = (1-w)*T_old(x,y) + w*1/4*(T(x-1,y)+T_old(x+1,y)+T(x,y-1)+T_old(x,y+1)); % SOR formula
                    end
                end
                err = max(abs(T-T_old),[],'all');
                cnt3(t)=cnt3(t)+1;
            end
        
            % Calculate FD
            for x = 2:m
                for y = 2:n
                    FD(x,y) = (T(x-1,y)+T(x+1,y)+T(x,y-1)+T(x,y+1)-4*T(x,y))/(dx*dx);
                end
            end
            T_mid(t) = T(m/2+1,n/2+1);
    
            % Plot filled contours with color and contour lines with labels
            contourf(T);
            colorbar;
            colormap(jet);
            
            % Add contour lines with labels
            hold on;
            contour(T, 'LineColor', 'k', 'ShowText', 'on');
            
            % Customize the plot
            xlabel('x');
            ylabel('y');
            title('SOR Solver');
    
            % Capture the current frame
            frame = getframe(gcf);
    
            % Convert the frame to an indexed image
            im = frame2im(frame);
            
            % Convert the indexed image to a colormap
            [imind, cm] = rgb2ind(im, 256);

            fileName = sprintf('outputs/%d/SORSolver.gif', i);
    
            imwrite(imind, cm, fileName, 'gif', 'WriteMode', 'append');
        
            % Pause between frames (optional)
            % pause(0.01);
    
            % Clear the figure for the next frame
            clf;
        end
    
        T_mid_3 = T_mid;
    
        % Plot of center point
        figure;
        plot(linspace(0,TT,timeStep+1),T_mid_1, 'b-', 'DisplayName', 'Jacobi Solver');
        hold on;
        
        plot(linspace(0,TT,timeStep+1),T_mid_2, 'g-', 'DisplayName', 'Gauss-Seidel Solver');
        hold on;
        
        plot(linspace(0,TT,timeStep+1),T_mid_3, 'r-', 'DisplayName', 'SOR Solver');
        
        % Customize the plot
        xlabel('t');
        ylabel('T_{mid}');
        title('Compare center point of each solver');
        
        % Add legend
        legend('show');

        fileName = sprintf('outputs/%d/CenterPointPlot.png', i);
    
        % Save the plot as a PNG file
        saveas(gcf, fileName);
    
        % Plot of iteration times
        figure;
        plot(linspace(0,TT,timeStep+1),cnt1, 'b-', 'DisplayName', 'Jacobi Solver');
        hold on;
        
        plot(linspace(0,TT,timeStep+1),cnt2, 'g-', 'DisplayName', 'Gauss-Seidel Solver');
        hold on;
        
        plot(linspace(0,TT,timeStep+1),cnt3, 'r-', 'DisplayName', 'SOR Solver');
        
        % Customize the plot
        xlabel('t');
        ylabel('Number of iteration times');
        title('Compare iteration times of each solver');
        
        % Add legend
        legend('show');

        fileName = sprintf('outputs/%d/IteNumPlot.png', i);
    
        % Save the plot as a PNG file
        saveas(gcf, fileName);
        clf;
        % fprintf('Testcase %d:\n', i);
        % fprintf('%f %f %f %f %f %f\n', T_mid_1(1),T_mid_2(1),T_mid_3(1),T_mid_1(timeStep+1),T_mid_2(timeStep+1),T_mid_3(timeStep+1));
        % fprintf('%d %d %d\n',cnt1(timeStep+1),cnt2(timeStep+1),cnt3(timeStep+1));
        clear;
    end
end