n_rocas=input("ingrese el numero de rocas: ");
for r=1:n_rocas
    clf, clearvars, clear
    g = -9.81;
    Vo = input("ingrese la velocidad inicial en m/s: "); 
    dt = 0.2; 

    %investigacion
    cd = 1.01; %coeficiente de drag
    densidad_del_aire = 1.29; 
    densidad_de_roca = 2600; 
    rad_roca = input("cual es el radio de la roca en m: ");
    flujo_del_aire = 3; 

    % Calculos
    area_de_roca = pi*rad_roca^2;
    volumen_de_roca = (4/3)*pi*rad_roca^3;
    masa_de_roca = densidad_de_roca*volumen_de_roca;
    b = (cd*densidad_del_aire*flujo_del_aire^2*area_de_roca)/(2*masa_de_roca); %constante de drag

    for a=20:5:50

        angulo = a;
        x = []; % vector de pociciones en x 
        y = []; % vector de posiciones en y 
        Vx = []; % velocidad en x con respecto al tiempo
        Vy = []; % velocidad en y con respecto al tiempo
        ax = [];
        ay = [];
        Vo_x = Vo.*cosd(angulo); %velocidad en x
        Vo_y = Vo.*sind(angulo); % velocidad inicial en y
        Vx(1) = Vo_x;
        Vy(1)= Vo_y;
        y(1)=463;
        x(1)= 0;
        ax(1) = 0;
        ay(1) = g;
        xn = x(1);  
        yn = y(1);
        xnm1 = xn - Vo*cosd(angulo)*dt;
        ynm1 = yn - Vo*sind(angulo)*dt + 0.5*g*dt^2;
        i = 1;

        while y(i) > 0
        % metodo de verlet
     
        xnp1 = (2*xn-xnm1)-.5*((xn-xnm1)/dt)^2*b/masa_de_roca*dt^2;
        x(i+1) = xnp1;
        ynp1 = (2*yn-ynm1) - (((yn-ynm1)/dt)/(abs((yn-ynm1)/dt)))*((yn-ynm1)/dt)^2*b/masa_de_roca*dt^2+.5*g*dt^2;
        y(i+1) = ynp1;
        Vx(i+1) = (x(i+1)-x(i))/dt;
        Vy(i+1) = (y(i+1)-y(i))/dt;
        ax(i+1) = (-b.*(Vx(i))^2) / masa_de_roca;
        ay(i+1) = (Vy(i)/abs(Vy(i)))* ((-b.*(Vy(i))^2) / masa_de_roca) + g;
        ynm1 = yn;
        xnm1 = xn;
        xn = xnp1;
        yn = ynp1;

        i = i+1;
    
        end
        
        xlim([0 5000]);
        ylim([0 2000]);

        yline(463,'--','altura del volcan')
        yline(0,'--','base del volcan')

        title('piroclasticos lanzados(20° a 50°)')
        xlabel('x')
        ylabel('y')

        hold on;
        xd = plot(x(1:i),y(1:i),'--');
        for i = 1:length(x)
            xd.XData = x(1:i);
            xd.YData = y(1:i);
            drawnow
            pause(0.0001)
        end
    end
end