%Limpieza de pantalla
clear all
close all
clc

%Declaración de variables simbólicas
syms th1(t) l1 

%Configuración del robot, 0 para junta rotacional, 1 para junta prismática
RP=[0];


%Creamos el vector de coordenadas generalizadas
Q= [th1];
disp('Coordenadas generalizadas');
pretty (Q);

%Creamos el vector de velocidades generalizadas
Qp= diff(Q, t);
disp('Velocidades generalizadas');
pretty (Qp);

%Número de grado de libertad del robot
GDL= size(RP,1);


%Vector de translación de la articulación 1 respecto a 0
P(:,:,1)= [l1*cos(th1); 
           l1*sin(th1);
                    0];

%Matriz de rotación de la junta 1 respecto a 0
R(:,:,1)= [cos(th1) -sin(th1)  0;
           sin(th1)  cos(th1)  0;
           0         0         1];


%Derivadas parciales de x respecto a th1
Jv11= functionalDerivative(P(1,1,GDL), th1);

%Derivadas parciales de y respecto a th1 
Jv21= functionalDerivative(P(2,1,GDL), th1);

%Derivadas parciales de z respecto a th1 
Jv31= functionalDerivative(P(3,1,GDL), th1);

disp('Cinemática diferencial de la posición del péndulo');

%Obtenemos la cinemática diferencial del péndulo a partir de la cinemática directa
jv_d=simplify([Jv11 ;
               Jv21 ; 
               Jv31]);

pretty(jv_d)

