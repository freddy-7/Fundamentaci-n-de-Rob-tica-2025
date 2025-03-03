%Limpieza de pantalla
clear all
close all
clc

%Declaración de variables simbólicas
syms th1(t) th2(t) th3(t) th4(t) th5(t) th6(t) t l1 l2 l3 l4 l5 l6

%Configuración del robot, 0 para junta rotacional, 1 para junta prismática
RP=[0 0 0 0 0 0];

%Creamos el vector de coordenadas articulares
Q= [th1, th2, th3, th4, th5, th6];
%disp('Coordenadas generalizadas');
%pretty (Q);

%Creamos el vector de velocidades generalizadas
Qp= diff(Q, t);
%disp('Velocidades generalizadas');
%pretty (Qp);
%Número de grado de libertad del robot
GDL= size(RP,2);
GDL_str= num2str(GDL);

%  rotacion_z= [cos(th1)  -sin(th1)     0 ;
%               sin(th1)   cos(th1)     0 ;
%               0             0         1];
%
%  rotacion_y= [cos(th1)    0      sin(th1) ;
%                  0        1          0    ;
%              -sin(th1)    0      cos(th1)];
%
%  rotacion_x= [1           0          0   ;
%               0        cos(th1)   -sin(th1);
%               0        sin(th1)   cos(th1)];
%
%  x_transf= [1  0  0;   %x +90   
%             0  0 -1;
%             0  1  0];
%            
%            [1  0  0;   %x -90
%             0  0  1;
%             0 -1  0];

%  y_transf= [0  0  1;   %y +90
%             0  1  0;
%            -1  0  0]; 


%Articulación 1
%Articulación 1 a Articulación 2
%Posición de la articulación 1 a 2
P(:,:,1)= [0;0;l1];
%Matriz de rotación de la junta 1 a 2---Transformación= Rot z(th1)*Rot y(+90)
R(:,:,1)= [0 -sin(th1) cos(th1);
           0  cos(th1) sin(th1);
          -1      0        0  ];

%Articulación 2
%Articulación 2 a Articulación 3
%Posición de la articulación 2 a 3
P(:,:,2)= [-l2*sin(th2); l2*cos(th2);0];
%Matriz de rotación de la junta 2 a 3---Transformación= Rot z(th2)
R(:,:,2)= [cos(th2) -sin(th2)  0;
           sin(th2)  cos(th2)  0;
           0         0         1];

%Articulación 3 
%Articulación 3 a Articulación 4
%Posición de la articulación 3 a 4
P(:,:,3)= [-l3*sin(th3); l3*cos(th3);0];
%Matriz de rotación de la junta 3 a 4---Transformación= Rot z(th3)*Rot x(-90)
R(:,:,3)=   [cos(th3)   0   -sin(th3);
             sin(th3)   0    cos(th3);
              0        -1        0  ];

%Articulación 4 
%Articulación 4 a Articulación 5
%Posición de la articulación 4 a 5
P(:,:,4)= [0; 0; l4];
%Matriz de rotación de la junta 4 a 5---Transformación= Rot z(th4)*Rot x(+90)
R(:,:,4)=   [cos(th4)   0    sin(th4);
             sin(th4)   0   -cos(th4);
              0         1        0  ];

%Articulación 5 
%Articulación 5 a Articulación 6
%Posición de la articulación 5 a 6
P(:,:,5)= [-l5*sin(th5); l5*cos(th5);0];
%Matriz de rotación de la junta 5 a 6---Transformación= Rot z(th5)*Rot y(+90)
R(:,:,5)=   [ 0    sin(th5)   cos(th5);
              0    cos(th5)   sin(th5);
             -1       0          0   ];

%Articulación 6 
%Articulación 6 a Extremo Final
%Posición de la articulación 6 a Extremo Final
P(:,:,6)= [0; 0; l6];
%Matriz de rotación de la junta 5 a 6---Transformación= Rot z(th6)
R(:,:,6)= [cos(th6) -sin(th6)  0;
           sin(th6)  cos(th6)  0;
           0         0         1];


%Creamos un vector de ceros
Vector_Zeros= zeros(1, 3);

%Inicializamos las matrices de transformación Homogénea locales
A(:,:,GDL)=simplify([R(:,:,GDL) P(:,:,GDL); Vector_Zeros 1]);
%Inicializamos las matrices de transformación Homogénea globales
T(:,:,GDL)=simplify([R(:,:,GDL) P(:,:,GDL); Vector_Zeros 1]);
%Inicializamos las posiciones vistas desde el marco de referencia inercial
PO(:,:,GDL)= P(:,:,GDL); 
%Inicializamos las matrices de rotación vistas desde el marco de referencia inercial
RO(:,:,GDL)= R(:,:,GDL); 


for i = 1:GDL
    i_str= num2str(i);
   %disp(strcat('Matriz de Transformación local A', i_str));
    A(:,:,i)=simplify([R(:,:,i) P(:,:,i); Vector_Zeros 1]);
   pretty (A(:,:,i))

   %Globales
    try
       T(:,:,i)= T(:,:,i-1)*A(:,:,i);
    catch
       T(:,:,i)= A(:,:,i);
    end
    disp(strcat('Matriz de Transformación global T', i_str))
    T(:,:,i)= simplify(T(:,:,i))
    pretty(T(:,:,i))

    RO(:,:,i)= T(1:3,1:3,i);
    PO(:,:,i)= T(1:3,4,i);
    %pretty(RO(:,:,i));
    %pretty(PO(:,:,i));
end


%Calculamos el jacobiano lineal de forma diferencial
%disp('Jacobiano lineal obtenido de forma diferencial');
%Derivadas parciales de x respecto a th1 y th2
Jv11= functionalDerivative(PO(1,1,GDL), th1);
Jv12= functionalDerivative(PO(1,1,GDL), th2);
Jv13= functionalDerivative(PO(1,1,GDL), th3);
%Derivadas parciales de y respecto a th1 y th2
Jv21= functionalDerivative(PO(2,1,GDL), th1);
Jv22= functionalDerivative(PO(2,1,GDL), th2);
Jv23= functionalDerivative(PO(2,1,GDL), th3);
%Derivadas parciales de z respecto a th1 y th2
Jv31= functionalDerivative(PO(3,1,GDL), th1);
Jv32= functionalDerivative(PO(3,1,GDL), th2);
Jv33= functionalDerivative(PO(3,1,GDL), th3);

%Creamos la matríz del Jacobiano lineal
jv_d=simplify([Jv11 Jv12 Jv13;
              Jv21 Jv22 Jv23;
              Jv31 Jv32 Jv33]);
%pretty(jv_d);


%Calculamos el jacobiano lineal de forma analítica
Jv_a(:,GDL)=PO(:,:,GDL);
Jw_a(:,GDL)=PO(:,:,GDL);

for k= 1:GDL
    if RP(k)==0 %Casos: articulación rotacional
       %Para las juntas de revolución
        try
            Jv_a(:,k)= cross(RO(:,3,k-1), PO(:,:,GDL)-PO(:,:,k-1));
            Jw_a(:,k)= RO(:,3,k-1);
        catch
            Jv_a(:,k)= cross([0,0,1], PO(:,:,GDL));%Matriz de rotación de 0 con respecto a 0 es la Matriz Identidad, la posición previa tambien será 0
            Jw_a(:,k)=[0,0,1];%Si no hay matriz de rotación previa se obtiene la Matriz identidad
         end
     %Para las juntas prismáticas
     elseif RP(k)==1 %Casos: articulación prismática
%         %Para las juntas prismáticas
        try
            Jv_a(:,k)= RO(:,3,k-1);
        catch
            Jv_a(:,k)=[0,0,1];
        end
            Jw_a(:,k)=[0,0,0];
     end
 end    

Jv_a= simplify (Jv_a);
Jw_a= simplify (Jw_a);
%disp('Jacobiano lineal obtenido de forma analítica');
%pretty (Jv_a);
%disp('Jacobiano ángular obtenido de forma analítica');
%pretty (Jw_a);


disp('Velocidad lineal obtenida mediante el Jacobiano lineal');
V=simplify (Jv_a*Qp');
pretty(V);
disp('Velocidad angular obtenida mediante el Jacobiano angular');
W=simplify (Jw_a*Qp');
    pretty(W);