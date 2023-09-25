clear all;
close all;
%Introducir datos
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\EXPERIMENTOS\SIN TIEMPO DE ESPERA\sensor_data_2023_3_3_12_9_300MBAR_10ms.csv');
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\EXPERIMENTOS\SIN TIEMPO DE ESPERA\sensor_data_2023_3_3_12_18_200MBAR_10ms.csv');
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\09_05\sensor_data_2023_4_2_10_44_Calibrado_SEspera_400mbar_55%.csv');
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\EXPERIMENTOS\SIN TIEMPO DE ESPERA\400mbar\sensor_data_2023_4_2_11_17_400MBAR_45%.csv');
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\EXPERIMENTOS\SIN TIEMPO DE ESPERA\400mbar\sensor_data_2023_4_2_11_33_400MBAR_50%.csv');
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\EXPERIMENTOS\SIN TIEMPO DE ESPERA\400mbar\sensor_data_2023_4_2_11_41_400MBAR_52.5%.csv');
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\EXPERIMENTOS\SIN TIEMPO DE ESPERA\400mbar\sensor_data_2023_3_4_13_52_400MBAR_10ms_55%.csv');
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\EXPERIMENTOS\SIN TIEMPO DE ESPERA\400mbar\sensor_data_2023_4_2_11_47_400MBAR_57.5%.csv');
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\EXPERIMENTOS\SIN TIEMPO DE ESPERA\400mbar\sensor_data_2023_4_2_11_19_400MBAR_60%.csv');
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\EXPERIMENTOS\SIN TIEMPO DE ESPERA\400mbar\sensor_data_2023_4_2_11_50_400MBAR_62.5%.csv');
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\EXPERIMENTOS\SIN TIEMPO DE ESPERA\400mbar\sensor_data_2023_4_2_11_39_400MBAR_65%.csv');

datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\EXPERIMENTOS\SIN TIEMPO DE ESPERA\400mbar\sensor_data_2023_3_4_13_52_400MBAR_10ms_55%.csv');

% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\09_05\sensor_data_2023_4_2_10_44_Calibrado_SEspera_400mbar_55%.csv');
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\23_3\bueno.csv');

% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\EXPERIMENTOS\400MBAR\SOLO RESP\sensor_data_2023_3_3_14_36_50%VALV.csv');
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\EXPERIMENTOS\400MBAR\INSP-EXH_2S_3S\sensor_data_2023_3_3_16_55_50%VALV.csv');
%VERIFICAR
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\Valvula\PRUEBAS PARA K VALVULA\400MBAR\sensor_data_2023_4_3_10_53_50%.csv');
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\Valvula\PRUEBAS PARA K VALVULA\400MBAR\sensor_data_2023_4_3_11_44_51%.csv');
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\Valvula\PRUEBAS PARA K VALVULA\400MBAR\sensor_data_2023_4_3_10_55_52.5%.csv');
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\Valvula\PRUEBAS PARA K VALVULA\400MBAR\sensor_data_2023_4_3_11_46_53%.csv');
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\Valvula\PRUEBAS PARA K VALVULA\400MBAR\sensor_data_2023_4_3_11_48_54%.csv');
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\Valvula\PRUEBAS PARA K VALVULA\400MBAR\sensor_data_2023_4_3_11_3_55%.csv');
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\Valvula\PRUEBAS PARA K VALVULA\400MBAR\sensor_data_2023_4_3_11_5_57.5%.csv');
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\Valvula\PRUEBAS PARA K VALVULA\400MBAR\sensor_data_2023_4_3_11_7_60%.csv');
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\Valvula\PRUEBAS PARA K VALVULA\400MBAR\sensor_data_2023_4_3_11_9_62.5%.csv');

%NUEVO GLOBO
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\Valvula\PRUEBAS CON GLOBO NUEVO\400MBAR\sensor_data_2023_4_3_11_56_50%.csv');
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\Valvula\PRUEBAS CON GLOBO NUEVO\400MBAR\sensor_data_2023_4_3_12_0_51%.csv');
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\Valvula\PRUEBAS CON GLOBO NUEVO\400MBAR\sensor_data_2023_4_3_12_4_52.5%.csv');
%datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\Valvula\PRUEBAS CON GLOBO NUEVO\400MBAR\sensor_data_2023_4_3_12_8_53%.csv');
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\Valvula\PRUEBAS CON GLOBO NUEVO\400MBAR\sensor_data_2023_4_3_12_10_54%.csv');
% datos=readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\Valvula\PRUEBAS CON GLOBO NUEVO\400MBAR\sensor_data_2023_4_3_12_13_55%.csv');


% datos= readmatrix('C:\Users\lucas\Desktop\INGENIERIA INDUSTRIAL\TFG\Laboratorio\ventilator-com-server-main\src\csv_logs\CASCADA\sensor_data_2023_4_3_11_59_EXP1.csv');

subtitulo={'Factor de olvido','a1','a2','b1','b2','Tiempo de muestreo'};
def={'1','0','0','0','0','0.02'};
val=inputdlg(subtitulo,'Valores iniciales',1,def);
fo=str2double(val(1));
%Hemos puesto un valor predeterminado del factor de olvido de 1. El
%algorítmo RLS es sensible al factor de olvido, por lo que es importante
%ajustarlo correctamente. Cuanto mayor sea, mas valor le dará a los datos
%recientes y cuanto menor sea mas valor a datos antiguos. Normalmente se 
%elige un valor entre 0.95 y 0.99. Elegiremos 1, para no tener factor de
%olvido ya que el sistema hace cambio repentinos que queremos seguir.
a1_ini=str2double(val(2));
a2_ini=str2double(val(3));
b1_ini=str2double(val(4));
b2_ini=str2double(val(5));
T=str2double(val(6));
%Orden del sistema a a proximar, orden 2.
orden=2;
%Inicializamos theta
theta=[a1_ini,a2_ini,b1_ini,b2_ini]';

%Le restamos a todas las filas del vector de tiempos el valor del primer término para que empiece
%en 0
datos(:,7)=datos(:,7)-datos(1,7);

%Nombramos un vector de tiempos que será el normalizado a cada paso de
%0.01s (tiempo de muestreo). Tenemos que hacer esto porque el plc nos envia
%datos en diferentes tiempos no constantes.
tiempos=zeros(length(datos),1);
for k=1:length(datos)-1
    tiempos(k+1,1)=tiempos(k,1)+(datos(k+1,7)-datos(k,7));
end
%Lo expresamos en segundos.
tiempos=tiempos/1000;

%La entrada es el flujo de aire en (L/min) y la salida la presión en mbar.
%Interpolamos para obtener los valores a cada 0.01s (Tiempo de muestreo).
salida1=interp1(tiempos',datos(:,1)',[0:T:tiempos(end)]);
salida2=interp1(tiempos',datos(:,2)',[0:T:tiempos(end)]);
entrada=interp1(tiempos',datos(:,3)',[0:T:tiempos(end)]);
valvula=interp1(tiempos',datos(:,5)',[0:T:tiempos(end)]);
%Todo al sistema internacional
% salida1=salida1.*100;       %Pascales
% salida2=salida2.*100;
salida1=salida1;       %Pascales
salida2=salida2;
ent=entrada./60;            %L/s
% ent=entrada;
val=valvula./100;
% valvula=valvula./100;       %Tanto por ciento
valvula=valvula
t=[0:T:tiempos(end)];

%Cuando no hay flujo el sensor de flujo mide a veces valores negativos, eso
%nos puede dar problemas a la hora de aproximar el módelo. Ante una entrada
%negativa la aproximación decrece como si se tratará de una presión
%negativa. Con este problema genera una "desventaja negativa" que cuando el
%sistema comienza a aproximarse crea un offset respecto la salida real.
ent_2=zeros(length(ent),1);
for k=1:length(ent)
    if ent(k)<0
        ent_2(k)=0;
    else
        ent_2(k)=ent(k);
    end
end
ent=ent_2;
% ent=valvula;
%Usamos la media aritmética de la presión a la entrada y la salida del
%pulmón.
p_average=(salida1+salida2)/2;
p_average2 = zeros(length(p_average), 1);
%La salida debe de empezar en 0 para que el método de mínimos
%cuadrados recurrentes o recursivos funcione correctamente.
for k=1:length(p_average)-1
    p_average2(k+1)=p_average(k)-p_average(1,1);
end
sal=p_average2;
%Definimos inicialmente la matriz Pn. Los coeficientes de la diagonal
%principal deben de ser >>1
Pn_ini=10000;
%La matriz de covarianza es el doble del orden del sistema, en este caso 4x4
Pn=Pn_ini*eye(2*orden);
%Definimos un vector para su evolución de los parámetros
evol=theta;
%Otro vector para la evolución de su error, inicializandolo en 0
e=0;
em=0;

%Vemos el tamaño del archivo de datos cargado anteriormente
tam_matriz=size(datos);
tam=tam_matriz(1,1);
%En este caso al interpolar, el vector va a ser mayor que lo que hemos
%cargado en los datos.
tam=size(sal);
%Condiciones iniciales de entrada y salida nulas
m_ent=0;
m_sal=0;
%Filtrar entrada y salida
s1=0; s2=0; alpha=1; salfil=zeros(length(sal),1); entfil=zeros(length(salfil),1);
for k=1:tam
    s1=(alpha*sal(k))+(1-alpha)*s1;
    salfil(k)=s1;
    s2=(alpha*ent(k))+(1-alpha)*s2;
    if s2<0
        entfil(k)=0;
    else
    entfil(k)=s2;
    end
end
sal=salfil;
ent=entfil;
% ent=valvula;
%% Algoritmo de mínimos cuadrados recursivos (RLS) con factor de olvido

for i=2*orden:tam
    %En este bucle vamos actualizando los valores de las entradas y salidas
    %reales.
    for j=1:orden
        m_sal(j)=[-sal(i-j)];
        m_ent(j)=[ent(i-j)];
    end
    %Vector de entradas y salidas
    M=[m_sal m_ent]';
    %Matriz de corrección
    Ln=(Pn*M)/(fo+(M'*Pn*M));
    %actualizamos el error
    e=sal(i)-M'*theta;
    %Calculo de theta+1
    theta=theta+Ln*e
    %Actualización de matriz de covarianza
    Pn=(Pn - Ln*M'*Pn)/fo;
    %Vector que acumula el error en cada iteración
    em=[em,e];
    %Vector que almacena la evolución de los parámetros en cada iteración
    evol=[evol,theta];
end

%Identificamos el sistema
Gd=tf([theta(3) theta(4)],[1 theta(1) theta(2)],T)

Gc=d2c(Gd,'zoh')
%% Linealización de la válvula
syms rho phi Pin Pa0  Pv0  real
D=Pv0*phi;
Q=Pv0*phi*sqrt(2*(Pin-Pa0)/rho); %Ecuación del flujo por la válvula
a1=Q;
a2=diff(Q,Pa0);
a3=diff(Q,Pv0);
Pv0 = subs(0.72);
Pa0 = subs(1000);
rho = subs(1.2);
phi = subs(0.005);
Pin= subs(40000);
a1=double(eval(a1));
a2=double(eval(a2));
a3=double(eval(a3));
flujo=zeros(length(ent),1);
Pa=sal;
% valvula1=0.175*x+0.45;
for k=1:length(flujo)
    valvula(k)=valvula(k)/0.625;
    if valvula(k)>1
        valvula(k)=1;
    end
    if valvula(k)<=0
        flujo(k)=0;
    else
%        flujo(k)=Pv0*phi*sqrt(-(2*Pa0-2*Pin)/rho)-((Pv0*phi)/(rho*sqrt(-(2*Pa0-2*Pin)/rho))*(Pa(k)-Pa0))+phi*sqrt(-(2*Pa0-2*Pin)/rho)*(valvula(k)-Pv0);
       flujo(k)=a1+a2*(Pa(k)-Pa0)+a3*(valvula(k)-Pv0);
    end
end
%% Gráficas de los resultados
[sal_sim,t_sim]=lsim(Gd,ent); %Entrada real de flujo para el sistema aproximado RLS
sal_simPOR=lsim(Gd,flujo);    %Entrada aproximada mediante la linealización de la válvulo para el sistema aproximado RLS
sal_simPOR=lsim(Gd,valvula); 
figure;
subplot(2,2,1:2)
plot(t,sal,t,sal_sim)
xlim([3 t(end)])
hold on;
plot(t,0.65*sal_simPOR,'Color',[0.4660 0.6740 0.1880]) 
xlim([3 t(end)])
hold off;
xlabel('Tiempo [s]') 
ylabel('Presion [mbar]')
legend('Sistema Real','Sistema estimado','Sistema estimado tomando la válvula como ganancia proporcional','Location','best','Orientation','horizontal')
subplot(2,2,3)
plot(t_sim,ent);
title('Entrada flujo ');
ylabel('Flujo [L/s]')
xlabel('Tiempo [s]');
xlim([3 t(end)])
subplot(2,2,4)
plot(t_sim,val,'Color',[0.4660 0.6740 0.1880]);
xlim([3 t(end)])
title('Entrada válvula ');
ylabel('Apertura válvula [%]')
xlabel('Tiempo [s]');



% figure;
% subplot(2,2,1);
% plot(t,sal,'b')
% legend('Señal real')
% title('Salida real (azul)');
% ylabel('Presión [mbar]');
% xlabel('Tiempo [seg.]');
% subplot(2,2,2);
% plot(t,sal_sim,'r');
% title('Salida del sistema estimado (rojo)');
% ylabel('Presión [mbar]');
% xlabel('Tiempo [seg.]');
% subplot(2,2,1:2);
% plot(t,sal,t,sal_sim,'r','LineWidth',1);
% title('Comparación entre la salida real y la salida del sistema estimado');
% legend('Salida Real','Salida Simulada','Location','best','Orientation','horizontal')
% ylabel('Presión [mbar]');
% xlabel('Tiempo [s]');
% xlim([2 t(end)])
% subplot(2,2,3:4);
% plot(t_sim,ent);
% title('Entrada real ');
% ylabel('Flujo [L/s]')
% xlabel('Tiempo [s]');
% xlim([2 t(end)])


% subplot(3,2,5:6)
% plot(t,[0 0 em])
% title('error cometido para cada iteración');
% xlabel('Tiempo [s]');
% ylabel('Error cuadrático medio')
% xlim([3 t(end)])
% subplot(3,2,5:6)
% plot(t,val)
% title('Porcentaje de apertura de la válvula')
% ylabel('Apertura [%]')
% xlabel('Tiempo [s]')
% xlim([2 t(end)])
% figure()
% plot(evol(1,:))
% hold on;
% plot(evol(2,:))
% hold on;
% plot(evol(3,:))
% hold on;
% plot(evol(4,:))
% hold off;
% title('Evolución de los parámetros estimados');
% xlabel('Muestras');
% ylabel('Evolución de los parámetros estimados')
% legend('a1','a2','b1','b2')
% 
% plot(t,valvula)
% title('Porcentaje de apertura de válvula')
% ylabel('Porcentaje %')
% xlabel('Tiempo [seg.]')
% figure()
% plot(t,ent,t,flujo)
% title('Comparación flujo medido y flujo calculado')
% legend({'Flujo real','Flujo aproximado'},'Location','Best','Orientation','horizontal')
% %Quedaría despejar los parámetros necesarios a partir de las formulas ya
% %conocidas de las matrices, o bien, directamente desde la función de
% %transferencia

% syms Em Rm E T z t C1 C2 tau real;
% A=[-Em/Rm Em*E/Rm; 0 0]; B=[Em+E;1]; C=[1 0]; D=0;
% H=int(expm(A*(t-tau))*B,tau,0,t);
% G=expm(A*t);
% H=subs(H,[E Em],[1/C1 1/C2]);
% G=subs(G,[E Em],[1/C1 1/C2]);
% Gdt=C*inv(z*[1 0;0 1]-G)*H;
% simplify(Gdt);
% coeficientes=collect(Gdt,z);
% eq1= (-C1-C2*exp(T/(C2*Rm)))/(C1*exp(T/(C2*Rm))) == theta(1);
% eq2= C1/(C1*exp(T/(C2*Rm)))== theta(2);
% eq3= (T*exp(T/(C2*Rm))-C1*Rm+C1*Rm*exp(T/(C2*Rm)))/(C1*exp(T/(C2*Rm)))== theta(3);
% eq4= (C1*Rm-T-C1*Rm*exp(T/(C2*Rm)))/(C1*exp(T/(C2*Rm)))== theta(4);
% % ecuaciones= subs([eq1 eq2 eq3 eq4],[C2 Rm T],[0.025 20 0.01])
% C2=solve(eq2,C2)
% C2=double(C2)
% % %Solo es necesaria la primera ecuación
% % eq=subs(eq1,[C2 Rm T],[0.025 20 0.05]);
% % C1=double(solve(eq,C1))
% % C1=strcat(num2str(C1),' L/mbar')