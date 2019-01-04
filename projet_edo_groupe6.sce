printf('\n\n                   *********** Soyez les Bienvenus************                         \n');
printf('\n       Projet-EDO-Groupe6-KOUAGOU N. Jean et BALOGOUN O.A. Ismaïla         \n\n');
                 printf('\n\n                  Début                  \n\n');

printf('\n\nToutes les exercices sont traités, veuillez consulter les figures et le fichier pdf \n\n');

printf('\n\n Les figures du numéro 0 au numéro 4 concernent le premier exercice, du numéro 5 au numéro 7 pour l''exercice 2, du numéro 8 au numéro 12 pour l''exercice 3, du numéro 13 au numéro 25 pour l''exercice 4 et le reste pour l''exercice 5\n\n')

//1. Équation logistique u'(t)=au(t)(1-bu(t))
// Schéma d'Euler explicite pour cette équation
// Définition de la fontion f(t,u(t))
function y=f(t,u)
    y=a*u.*(1-b*u);
endfunction
//Définition de la fonction solEuler qui calcule la solution avec le schéma d'Euler explicite
function [solutionE]=solEulerE(temps,u0,f)
    solutionE=u0;
     h=temps(2)-temps(1);
    u=u0;
    for inter=1:length(temps)-1
        u=u+h*f(inter,u);
        solutionE=[solutionE,u]
    end
endfunction

//Définition de la fonction solEulerM qui calcule la solution avec le schéma d'Euler modifié
function [solutionEm]=solEulerM(temps,u0,f)
     h=temps(2)-temps(1);
    solutionEm=u0;
    u=u0;
    for inter=1:length(temps)-1
       ff=h/2*f(inter,u);
        u=u+h*f(inter+h/2,u+ff);
        solutionEm=[solutionEm,u]
    end
endfunction

// 2. Système de Lotka-Volterra(proies et prédateurs)
//Définition de la fonction g(t,y(t))=(ay1-by1y2,-cy2+dy1y2) avec y1 et y2 les composantes de y 
function u=g(t,y)
    u(1)=a*y(1)-b*y(1)*y(2);
    u(2)=-c*y(2)+d*y(1)*y(2);
endfunction

//Définition de la fonction solEuler1 qui calcule la solution avec le schéma d'Euler explicite
function [solutionE]=solEuler1(temps,y0,g)
    solutionE=y0;
     h=temps(2)-temps(1);
    y=y0;
    for inter=1:length(temps)-1
        y=y+h*g(inter,y);
        solutionE=[solutionE,y]
    end
endfunction

//Définition de la fonction solRK4 qui calcule la solution avec le schéma de RK4
function [solution]=solRK4(temps,y0,g)
    h=temps(2)-temps(1);
    solution=y0;
    s=y0;
    for t=temps(1:$-1)
    yn1=s;
    yn2=yn1+0.5*h*g(t,yn1);
    yn3=yn1+0.5*h*g(t+0.5*h,yn2);
    yn4=yn1+h*g(t+0.5*h,yn3);
    s=yn1+h/6*(g(t,yn1)+2*g(t+0.5*h,yn2)+2*g(t+0.5*h,yn3)+g(t+h,yn4));
    solution=[solution,s];
   
    
 end
endfunction

//Définition de la fonction CrankNicholson qui calcule la solution avec la méthode de Crank-Nicholson (RK implicite choisi)
//définition de la fonction reso qui sera utilisée par la fonction CrankNicholson 
function y=reso(x)
    y=x-x0-h/2*(g(t,x0)+g(t+h,x));
endfunction

function [solutionC]=CrankNicholson(temps,y0,g)
    solutionC=y0;
     h=temps(2)-temps(1);
    x0=y0;
    for t=1:length(temps)-1
        yy=fsolve(x0,reso);
        solutionC=[solutionC,yy];
        x0=yy;
    end
endfunction

//3-Modèle logistique de Verhuls
//définition de la fonction associée au système
function u=k(t,y)
    u(1)=a*y(1)-e*(y(1))^2-b*y(1)*y(2);
    u(2)=-c*y(2)+d*y(1)*y(2);
endfunction

//4-Populatios en compétition
//définition de la fonction associé au système 
function u=j(t,y)
    u(1)=y(1)*(1-y(1)-(o12)*y(2));
     u(2)=a*y(2)*(1-y(2)-(o21)*y(1));
endfunction

//// definition des fonctions seulement utilisables a la consigne 4 Compétition de deux espèces. il s'agit des droites de l'équilibre intérieur 

function y=f1(t)
    y=-(1/o12)*t+1/o12
endfunction
function y=f2(t)
    y=-(o21)*t+1
endfunction

//// definition de la fonctions utilisable a la consigne 5 généralisation au cas de trois espèces
function u=G(t,y)
    u(1)=a*y(1)-b*y(1)*y(2)-c*y(1)*y(3);
    u(2)=-d*y(2)-e*y(2)*y(3)+ff*y(1)*y(2);
    u(3)=-gg*y(3)-hh*y(2)*y(3)+i*y(1)*y(3);
endfunction
//programme principal pour Euler explicite 
//(a,b)=(0.5,0.1)
a=0.5;
b=0.1;
T=1;
temps=linspace(0,T);

xset('window',0)
//Affichage des solutions pour différentes valeurs de u0 avec (a,b)=(0.5,0.1)
clf();
for u0=0:10
    [solutionE]=solEulerE(temps,u0,f);
    plot(temps,solutionE)
end
title('Solution euler explicite pour (a,b)=(0.5,0.1) et u0 varie')

//(a,b)=(7,0.14)
a=7;
b=0.14;
//Affichage des solutions pour différentes valeurs de u0 avec (a,b)=(7,0.14)
xset('window',1)
clf();
for u0=0:10
    [solutionE]=solEulerE(temps,u0,f);
    plot(temps,solutionE)
end
xlabel('temps')
ylabel('évolution des populations')
xtitle('Solution euler explicite pour (a,b)=(7,0.14) et u0 varie')

//programme principal pour euler modifieé
//(a,b)=(0.5,0.1)
a=0.5;
b=0.1;
u0=4;
//Affichage de la solution avec Euler modifié pour a=0.5, b=0.1
xset('window',2)
clf
for u0=0:10
    [solutionEm]=solEulerM(temps,u0,f);
    plot(temps,solutionEm)
end
xlabel('temps')
ylabel('évolution des populations')
xtitle('Solution euler modifier pour (a,b)=(0.5,0.1) et u0 varie')


//(a,b)=(7,0.14)
a=7;
b=0.14;
//Affichage des solutions avec Euler modifié pour a=7, b=0.14
xset('window',3)
clf
for u0=0:10
    [solutionEm]=solEulerM(temps,u0,f);
    plot(temps,solutionEm)
end
xlabel('temps')
ylabel('évolution des populations')
xtitle('Solution euler modifié pour (a,b)=(7,0.14) et u0 varie')
//solution donnée par le solveur de scilab avec 
//(a,b)=(0.5,0.1)
a=0.5;
b=0.1;
u0=4;
[solode]=ode(u0,0,temps,f)
//comparaison de ode avec Euler explicite dans la même fenêtre

xset('window',0)

plot(temps,solode,'r+')

//comparaison de ode avec Euler modifié dans la même fenêtre

xset('window',2)

plot(temps,solode,'r+')
//Autre valeur pou a et b
//(a,b)=(7,0.14)
a=7;
b=0.14;
u0=4;

[solode]=ode(u0,0,temps,f)

//comparaison avec Euler Explicite

xset('window',1)

plot(temps,solode,'r+')

//comparaison avec Euler modifié

xset('window',3)

plot(temps,solode,'r+')

// on fait varier a

xset('window',4)

clf();
b=0.14;
subplot(1,2,1)
u0=1;
for a=0.5:0.5:10
    [solutionEm]=solEulerM(temps,u0,f);
    plot(temps,solutionEm)
end
xtitle('Solution euler modifié pour u0=1 b=0.14 et a varie')
subplot(1,2,2)

u0=10;
for a=0.5:0.5:10
    [solutionEm]=solEulerM(temps,u0,f);
    plot(temps,solutionEm)
end
xtitle('Solution euler modifié pour u0=10 b=0.14 et a varie')
// 2.Système de Lotka-Volterra(proies et prédateurs)
// Programme principale Euler explicite
a=1;
b=0.2;
c=0.5;
d=1;
T=1;
temps=linspace(0,T);

y0=[50;20];

[solutionE]=solEuler1(temps,y0,g);

xset('window',5)

clf();
subplot(2,2,1)
plot(temps,solutionE(1,:),'r')
plot(temps,solutionE(2,:))

//y0=[0;20]

y0=[0;20];
[solutionE]=solEuler1(temps,y0,g);
subplot(2,2,2)
plot(temps,solutionE(1,:),'r')
plot(temps,solutionE(2,:))

//y0=[50;0]

y0=[50;0];
[solutionE]=solEuler1(temps,y0,g);
subplot(2,2,3)
plot(temps,solutionE(1,:),'r')
plot(temps,solutionE(2,:))

//y0=[20;50]

y0=[20;50];
[solutionE]=solEuler1(temps,y0,g);
subplot(2,2,4)
plot(temps,solutionE(1,:),'r')
plot(temps,solutionE(2,:))

// Programme principal RK4
y0=[50;20];

[solution]=solRK4(temps,y0,g);

xset('window',6)

clf();
subplot(2,2,1)
plot(temps,solution(1,:),'r')
plot(temps,solution(2,:))

//y0=[0;20]

y0=[0;20];
[solution]=solRK4(temps,y0,g);
subplot(2,2,2)
plot(temps,solution(1,:),'r')
plot(temps,solution(2,:))

//y0=[50;0]

y0=[50;0];
[solution]=solRK4(temps,y0,g);
subplot(2,2,3)
plot(temps,solution(1,:),'r')
plot(temps,solution(2,:))

//y0=[20;50]

y0=[20;50];
[solution]=solRK4(temps,y0,g);
subplot(2,2,4)
plot(temps,solution(1,:),'r')
plot(temps,solution(2,:))

// Programme principal Crank-Nicholson

y0=[50;20];

[solutionC]=CrankNicholson(temps,y0,g);

xset('window',7)

clf
subplot(2,2,1)
plot(temps,solutionC(1,:),'r')
plot(temps,solutionC(2,:))

//y0=[0;20]

y0=[0;20];
[solutionC]=CrankNicholson(temps,y0,g);
subplot(2,2,2)
plot(temps,solutionC(1,:),'r')
plot(temps,solutionC(2,:))

//y0=[50;0]

y0=[50;0];
[solutionC]=CrankNicholson(temps,y0,g);
subplot(2,2,3)
plot(temps,solutionC(1,:),'r')
plot(temps,solutionC(2,:))

//y0=[20;50]

y0=[20;50];
[solutionC]=CrankNicholson(temps,y0,g);
subplot(2,2,4)
plot(temps,solutionC(1,:),'r')
plot(temps,solutionC(2,:))

//solution donnée par le solveur de scilab avec y0=[50,20]

y0=[50;20];
[solode]=ode(y0,0,temps,g);

//comparaison avec Euler Explicite

xset('window',5)

subplot(2,2,1)
plot(temps,solode(1,:),'o')
plot(temps,solode(2,:),'o')

//comparaison avec RK4

xset('window',6)

subplot(2,2,1)
plot(temps,solode(1,:),'o')
plot(temps,solode(2,:),'o')

//comparaison avec Crank-Nicholson

xset('window',7)
subplot(2,2,1)
plot(temps,solode(1,:),'o')
plot(temps,solode(2,:),'o')

//3.Modèle logistique de Verhuls
//Programme principale

//Euler explicite

temps=linspace(0,T);
y0=[50;20];

xset('window',8)
clf();

for e=[0,0.2,0.5,1] 
    [solutionE]=solEuler1(temps,y0,k);
    plot(temps,solutionE(1,:),'r')
    plot(temps,solutionE(2,:))

end
xlabel('temps')
ylabel('évolution des populations')
xtitle('y0 fixé, e varie')

//RK4

y0=[95;20];
xset('window',9)
clf();

for e=[0,0.2,0.5,1] 
    
    [solution]=solRK4(temps,y0,k);
    plot(temps,solution(1,:),'r')
    plot(temps,solution(2,:))
   
end
xlabel('temps')
ylabel('évolution des populations')
xtitle('y0 fixé e varie')

y0=[10;0];

xset('window',10)
clf();
for e=[5,7.5,9,11,14]
    
    [solution]=solRK4(temps,y0,k);
    subplot(2,1,1)
    plot(temps,solution(1,:),'r')
    plot(temps,solution(2,:))
end
for e=[0.001,0.002,0.005,0.009] 
    [solution]=solRK4(temps,y0,k);
    subplot(2,1,2)
    plot(temps,solution(1,:),'r')
    plot(temps,solution(2,:))
   end
   
   xtitle('y0=(10,0) fixé, e prenant les valeurs:0.001,0.002,0.005,0.009')

 
//Crank-Nicholson
y0=[30;80];

xset('window',11)
clf
for e=[0,0.2,0.5,1] 
    [solutionC]=CrankNicholson(temps,y0,k);
    plot(temps,solutionC(1,:),'r')
    plot(temps,solutionC(2,:))
   
end
xlabel('temps')
ylabel('évolution des populations')
xtitle('y0 fixé, e varie')

//solution donnée par le solveur de scilab avec e=1
y0=[50;20];
e=1;
[solode]=ode(y0,0,temps,k);

//comparaison avec Euler Explicite
xset('window',8)
plot(temps,solode(1,:),'o')
plot(temps,solode(2,:),'o')
//comparaison avec RK4
y0=[95;20];
xset('window',9)
plot(temps,solode(1,:),'o')
plot(temps,solode(2,:),'o')
//comparaison avec Crank-Nicholson
y0=[30;80];
xset('window',11)
plot(temps,solode(1,:),'o')
plot(temps,solode(2,:),'o')

//Etudions l'évolution du système en partant de la donnée initiale (a/e, 0).

//RK4
xset('window',12)
clf
 for e=[0.02,0.03,0.5,0.6,4]
    y0=[a/e;0]
    [solution]=solRK4(temps,y0,k);
    plot(temps,solution(1,:),'ro')
    plot(temps,solution(2,:))

end
xlabel('temps')
ylabel('évolution des populations')
xtitle('y=[a/e,0] et e varie')
//
//4.Populations en compétition
//Programme principale avec a=1, o12=3/2 et o21=1/2;y0=[50;20]
//Euler explicite
y0=[50;20];

xset('window',13)
clf();
subplot(2,3,1)
a=1;
o12=3/2;
o21=1/2;
    [solutionE]=solEuler1(temps,y0,j);
    plot(temps,solutionE(1,:),'r')
    plot(temps,solutionE(2,:))
xlabel('temps')
ylabel('évolution des populations')

//RK4

xset('window',14)
clf();
subplot(2,3,1)
    [solution]=solRK4(temps,y0,j)
    plot(temps,solution(1,:),'r')
    plot(temps,solution(2,:))
xlabel('temps')
ylabel('évolution des populations')
//Crank-Nicholson

xset('window',15)
clf
  subplot(2,3,1)
    [solutionC]=CrankNicholson(temps,y0,j);
    plot(temps,solutionC(1,:),'r')
    plot(temps,solutionC(2,:))
xlabel('temps')
ylabel('évolution des populations')

//solution donnée par le solveur de scilab avec y0=[50,20]

[solode]=ode(y0,0,temps,j);

//comparaison avec Euler Explicite

xset('window',13)
subplot(2,3,1)
plot(temps,solode(1,:),'o')
plot(temps,solode(2,:),'o')
//comparaison avec RK4
xset('window',14)
subplot(2,3,1)
plot(temps,solode(1,:),'o')
plot(temps,solode(2,:),'o')

//comparaison avec Crank-Nicholson

xset('window',15)
subplot(2,3,1)
plot(temps,solode(1,:),'o')
plot(temps,solode(2,:),'o')


//Programme principale avec a=1, o12=3/2 et o21=1/2;y0=[50;50]
//Euler explicite

y0=[50;50];

xset('window',13)
subplot(2,3,2)
    [solutionE]=solEuler1(temps,y0,j);
    plot(temps,solutionE(1,:),'r')
    plot(temps,solutionE(2,:))


//RK4
xset('window',14)
subplot(2,3,2)
    [solution]=solRK4(temps,y0,j);
    plot(temps,solution(1,:),'r')
    plot(temps,solution(2,:))

//Crank-Nicholson
xset('window',15)
  subplot(2,3,2)
    [solutionC]=CrankNicholson(temps,y0,j);
    plot(temps,solutionC(1,:),'r')
    plot(temps,solutionC(2,:))
    
//Programme principale avec a=1, o12=3/2 et o21=1/2;y0=[50;70]
//Euler explicite
y0=[50;70];
xset('window',13)
subplot(2,3,3)
    [solutionE]=solEuler1(temps,y0,j);
    plot(temps,solutionE(1,:),'r')
    plot(temps,solutionE(2,:))


//RK4
xset('window',14)
subplot(2,3,3)
    [solution]=solRK4(temps,y0,j);
    plot(temps,solution(1,:),'r')
    plot(temps,solution(2,:))

//Crank-Nicholson
xset('window',15)
  subplot(2,3,3)
    [solutionC]=CrankNicholson(temps,y0,j);
    plot(temps,solutionC(1,:),'r')
    plot(temps,solutionC(2,:))
    
//Programme principale avec a=1, o12=3/2 et o21=1/2;y0=[0;20]
//Euler explicite
y0=[0;20];

xset('window',13)
subplot(2,3,4)
    [solutionE]=solEuler1(temps,y0,j);
    plot(temps,solutionE(1,:),'r')
    plot(temps,solutionE(2,:))


//RK4
xset('window',14)
subplot(2,3,4)
    [solution]=solRK4(temps,y0,j);
    plot(temps,solution(1,:),'r')
    plot(temps,solution(2,:))

//Crank-Nicholson
xset('window',15)
  subplot(2,3,4)
    [solutionC]=CrankNicholson(temps,y0,j);
    plot(temps,solutionC(1,:),'r')
    plot(temps,solutionC(2,:))

//Programme principale avec a=1, o12=3/2 et o21=1/2;y0=[20;0]
//Euler explicite
y0=[20;0];
xset('window',13)
subplot(2,3,5)
    [solutionE]=solEuler1(temps,y0,j);
    plot(temps,solutionE(1,:),'r')
    plot(temps,solutionE(2,:))


//RK4
xset('window',14)
subplot(2,3,5)
    [solution]=solRK4(temps,y0,j);
    plot(temps,solution(1,:),'r')
    plot(temps,solution(2,:))

//Crank-Nicholson
xset('window',15)
  subplot(2,3,5)
    [solutionC]=CrankNicholson(temps,y0,j);
    plot(temps,solutionC(1,:),'r')
    plot(temps,solutionC(2,:))
//
//Programme principale avec a=1/2, o12=2 et o21=3;y0=[50;20]
//Euler explicite
y0=[50;20];
xset('window',16)
clf
subplot(2,3,1)
a=1/2;
o12=2;
o21=3;
    [solutionE]=solEuler1(temps,y0,j);
    plot(temps,solutionE(1,:),'r')
    plot(temps,solutionE(2,:))

  xlabel('temps')
  ylabel('évolution des populations')
  
//RK4

xset('window',17)
clf
subplot(2,3,1)
    [solution]=solRK4(temps,y0,j);
    plot(temps,solution(1,:),'r')
    plot(temps,solution(2,:))
xlabel('temps')
ylabel('évolution des populations')

//Crank-Nicholson

xset('window',18)
clf();
  subplot(2,3,1)
    [solutionC]=CrankNicholson(temps,y0,j);
    plot(temps,solutionC(1,:),'r')
    plot(temps,solutionC(2,:))
xlabel('temps')
ylabel('évolution des populations')

//solution donnée par le solveur de scilab avec y0=[50,20]

[solode]=ode(y0,0,temps,j);
//comparaison avec Euler Explicite
xset('window',16)
subplot(2,3,1)
plot(temps,solode(1,:),'o')
plot(temps,solode(2,:),'o')

//comparaison avec RK4
xset('window',17)
subplot(2,3,1)
plot(temps,solode(1,:),'o')
plot(temps,solode(2,:),'o')

//comparaison avec Crank-Nicholson

xset('window',18)
subplot(2,3,1)
plot(temps,solode(1,:),'o')
plot(temps,solode(2,:),'o')

//Programme principale aveca=1/2, o12=2 et o21=3;y0=[50;50]

//Euler explicite

y0=[50;50];
xset('window',16)
subplot(2,3,2)
    [solutionE]=solEuler1(temps,y0,j);
    plot(temps,solutionE(1,:),'r')
    plot(temps,solutionE(2,:))


//RK4
xset('window',17)
subplot(2,3,2)
    [solution]=solRK4(temps,y0,j);
    plot(temps,solution(1,:),'r')
    plot(temps,solution(2,:))

//Crank-Nicholson
xset('window',18)
  subplot(2,3,2)
    [solutionC]=CrankNicholson(temps,y0,j);
    plot(temps,solutionC(1,:),'r')
    plot(temps,solutionC(2,:))
    
//Programme principale avec a=1/2, o12=2 et o21=3;y0=[50;70]
//Euler explicite

y0=[50;70];
xset('window',16)
subplot(2,3,3)
    [solutionE]=solEuler1(temps,y0,j);
    plot(temps,solutionE(1,:),'r')
    plot(temps,solutionE(2,:))


//RK4

xset('window',17)
subplot(2,3,3)
    [solution]=solRK4(temps,y0,j);
    plot(temps,solution(1,:),'r')
    plot(temps,solution(2,:))

//Crank-Nicholson

xset('window',18)
  subplot(2,3,3)
    [solutionC]=CrankNicholson(temps,y0,j);
    plot(temps,solutionC(1,:),'r')
    plot(temps,solutionC(2,:))
    
//Programme principale avec a=1/2, o12=2 et o21=3;y0=[0;20]

//Euler explicite

y0=[0;20];
xset('window',16)
subplot(2,3,4)
    [solutionE]=solEuler1(temps,y0,j);
    plot(temps,solutionE(1,:),'r')
    plot(temps,solutionE(2,:))


//RK4
xset('window',17)
subplot(2,3,4)
    [solution]=solRK4(temps,y0,j);
    plot(temps,solution(1,:),'r')
    plot(temps,solution(2,:))

//Crank-Nicholson
xset('window',18)
  subplot(2,3,4)
    [solutionC]=CrankNicholson(temps,y0,j);
    plot(temps,solutionC(1,:),'r')
    plot(temps,solutionC(2,:))

//Programme principale avec a=1/2, o12=2 et o21=3;y0=[20;0]
//Euler explicite
y0=[20;0];
xset('window',16)
subplot(2,3,5)
    [solutionE]=solEuler1(temps,y0,j);
    plot(temps,solutionE(1,:),'r')
    plot(temps,solutionE(2,:))


//RK4
xset('window',17)
subplot(2,3,5)
    [solution]=solRK4(temps,y0,j);
    plot(temps,solution(1,:),'r')
    plot(temps,solution(2,:))

//Crank-Nicholson

xset('window',18)

  subplot(2,3,5)
    [solutionC]=CrankNicholson(temps,y0,j);
    plot(temps,solutionC(1,:),'r')
    plot(temps,solutionC(2,:))
  //on change les valeurs de a, o12, o21  
    a=1/2;
o12=0;
o21=3;

//Euler explicite
y0=[20;50];
xset('window',19)
    [solutionE]=solEuler1(temps,y0,j);
    plot(temps,solutionE(1,:),'r')
    plot(temps,solutionE(2,:))
    xlabel('temps')
    ylabel('évolution des populations')
    xtitle('o12=0 et o21=3 : l`une seulement des espèces inhibe le développement de l`autre')

    

//RK4
xset('window',20)
    [solution]=solRK4(temps,y0,j);
    plot(temps,solution(1,:),'r')
    plot(temps,solution(2,:))
    xlabel('temps')
    ylabel('évolution des populations')

//Crank-Nicholson
xset('window',21)
    [solutionC]=CrankNicholson(temps,y0,j);
    plot(temps,solutionC(1,:),'r')
    plot(temps,solutionC(2,:))
    xlabel('temps')
    ylabel('évolution des populations')
    
    //on modifie les valeurs de a, o12, o21
    a=1.5;
o12=0;
o21=0;

//Euler explicite
y0=[20;50];

xset('window',22)

    [solutionE]=solEuler1(temps,y0,j);
    plot(temps,solutionE(1,:),'r')
    plot(temps,solutionE(2,:))

   xlabel('temps')
   ylabel('évolution des populations')
//RK4
xset('window',23)
    [solution]=solRK4(temps,y0,j);
    plot(temps,solution(1,:),'r')
    plot(temps,solution(2,:))
    
    xlabel('temps')
    ylabel('évolution des populations')
    
//Crank-Nicholson

xset('window',24)
    [solutionC]=CrankNicholson(temps,y0,j);
    plot(temps,solutionC(1,:),'r')
    plot(temps,solutionC(2,:))
    xlabel('temps')
    ylabel('évolution des populations')
    xtitle('o12=o21=0 :la présence de l`une des populations est équivalente à son absence pour l`autre')

// on modifie les valeurs de y0, T et h
y0=[30;50];
T=3;
h=0.005;
temps=0:h:T;

xset('window',25)
clf()

a=1;
o12=5/3;
o21=4;

    y0=[50;50]
 [solutionC]=CrankNicholson(temps,y0,j);
    plot(temps,solutionC(1,:),'r')
    plot(temps,solutionC(2,:))
    plot(temps,zeros(temps),'g')
    xlabel('temps')
    ylabel('évolution des populations')
    xtitle('Les deux droites se coupent, les deux espèces peuvent cohabiter')
    u1=f1(temps);
    plot(temps,u1,'r')
    
    u2=f2(temps);
    plot(temps,u2)
    
    
// 5.Généralisation

a=2;
b=0.2;
c=0.9;
d=0.5;
e=0.01;
ff=0.1;
gg=1.5;
hh=0.1;
i=0.6;
T=1.6;
temps=linspace(0,T);
y0=[500;15;10];
[solutionC]=CrankNicholson(temps,y0,G);
xset('window',26)
clf();
plot(temps,solutionC(1,:),'r')
plot(temps,solutionC(2,:),'g')
plot(temps,solutionC(3,:))
xtitle('u0(1)=500  ,u0(2)=15   et  u0(3)=10')
xlabel('temps')
ylabel('évolution des trois populations');

//on modifie y0, a, b, c, e, et h 
a=3;
b=0.5;
c=0.7;
d=0.5;
e=0.5;
ff=0.1;
gg=1.5;
hh=1;
i=0.6; 
y0=[155;5;12];
[solutionC]=CrankNicholson(temps,y0,G);
xset('window',27)
plot(temps,solutionC(1,:),'r')
plot(temps,solutionC(2,:),'g')
plot(temps,solutionC(3,:))
xlabel('temps')
ylabel('évolution des trois populations')
xtitle('u0(1)=155  u0(2)=5   et  u0(3)=12');
//on modifie y0, a, ff, i, et T
a=1;
b=0.5;
c=0.7;
d=0.5;
e=0.5;
ff=0.3;
gg=1.5;
hh=1;
i=0.7;
T=3;
temps=linspace(0,T);
y0=[0;16;25];
xset('window',28)
clf();
[solutionC]=CrankNicholson(temps,y0,G);
plot(temps,solutionC(1,:),'r+')
plot(temps,solutionC(2,:),'g')
plot(temps,solutionC(3,:),)
xlabel('temps')
ylabel('évolution des trois populations')
xtitle('u0(1)=0   u0(2)=16   et  u0(3)=25');

//on modifie y0, a, ff, T et i
a=3;
b=0.5;
c=0.7;
d=0.5;
e=0.5;
ff=0.1;
gg=1.5;
hh=1;
i=0.6;
T=1;
temps=linspace(0,T);
y0=[52;0;2];
xset('window',29)
clf();
[solutionC]=CrankNicholson(temps,y0,G);
plot(temps,solutionC(1,:),'r')
plot(temps,solutionC(2,:),'g+')
plot(temps,solutionC(3,:),)
xlabel('temps')
ylabel('évolution des trois populations')
xtitle('u0(1)=52   u0(2)=0   et  u0(3)=2');

//on modifie y0 et T
a=3;
b=0.5;
c=0.7;
d=0.5;
e=0.5;
ff=0.1;
gg=1.5;
hh=1;
i=0.6;
T=2;
temps=linspace(0,T);
y0=[55;3;0];
xset('window',30)
clf();
[solutionC]=CrankNicholson(temps,y0,G);
plot(temps,solutionC(1,:),'r')
plot(temps,solutionC(2,:),'g')
plot(temps,solutionC(3,:),'b+')
xlabel('temps')
ylabel('évolution des trois populations')
xtitle('u0(1)=55   u0(2)=3   et  u0(3)=0');
//on modifie y0 et T
a=3;
b=0.5;
c=0.7;
d=0.5;
e=0.5;
ff=0.1;
gg=1.5;
hh=1;
i=0.6;
T=5;
temps=linspace(0,T);

//une perturbation du point d'équilibre (g/i,0,a/c)
y0=[gg/i+1.5;2;a/c+1.5];
xset('window',31)
clf();
[solutionC]=CrankNicholson(temps,y0,G);
plot(temps,solutionC(1,:),'r')
plot(temps,solutionC(2,:),'g+')
plot(temps,solutionC(3,:),'b-')
xlabel('temps')
ylabel('évolution des trois populations')
xtitle('u0(1)=g/i +1.5   u0(2)=2   et  u0(3)=a/c +1.5');


//une perturbation du point d'équilibre (g/i,0,a/c)
y0=[d/ff+1.5;a/b+1.5;5];
xset('window',32)
clf();
[solutionC]=CrankNicholson(temps,y0,G);
plot(temps,solutionC(1,:),'r')
plot(temps,solutionC(2,:),'g+')
plot(temps,solutionC(3,:),'b-')
xlabel('temps')
ylabel('évolution des trois populations')
xtitle('u0(1)=d/f +1.5   u0(2)=a/b +1.5   et  u0(3)=5');

//stabilité du cinquième équilibre

a=3;
b=0.5;
c=0.7;
d=0.5;
e=0.5;
ff=0.1;
gg=1.5;
hh=1;
i=0.6;

// la matrice jacobienne A de la fonction G au cinquième point d'équilibre
A=[0,(a*hh*e*b-e*b*b*gg-hh*b*d*c)/(hh*c*ff+e*i*b),(-e*c*b*gg-hh*d*c*c-a*hh*c*e)/(hh*c*ff+e*i*b);(i*d*c*ff-ff*ff*gg*c+a*ff*e*i)/(hh*c*ff+e*i*b),0,(-e*i*d*c+e*ff*gg*c-a*e*e*i)/(hh*c*ff+e*i*b);(-i*i*d*b+ff*gg*b*i+a*hh*ff*i)/(hh*c*ff+e*i*b),(i*b*d*hh-ff*gg*b*hh-a*hh*hh*ff)/(hh*c*ff+e*i*b),0];

T=5;
temps=linspace(0,T);

//une perturbation du point d'équilibre (g/i,0,a/c)

y0=[(e*b*gg+hh*d*c+a*hh*e)/(hh*c*ff+e*i*b);(i*d*c-ff*gg*c+a*e*i)/(hh*c*ff+e*i*b);(ff*gg*b+a*hh*ff-i*d*b)/(hh*c*ff+e*i*b)];
//affichage des figures avec Cranck Nicolson
xset('window',33)
clf();
[solutionC]=CrankNicholson(temps,y0,G);
plot(temps,solutionC(1,:),'r')
plot(temps,solutionC(2,:),'g+')
plot(temps,solutionC(3,:),'b-')
xlabel('temps')
ylabel('évolution des trois populations')
xtitle('cinquième état d''équilibre');

//perturbation du cinquième état d''équilibre

y0=[(e*b*gg+hh*d*c+a*hh*e)/(hh*c*ff+e*i*b)-1;(i*d*c-ff*gg*c+a*e*i)/(hh*c*ff+e*i*b)+1;(ff*gg*b+a*hh*ff-i*d*b)/(hh*c*ff+e*i*b)+0.5];

//on modifie T
T=15;
temps=linspace(0,T);
//affichage des figures avec Cranck Nicolson
xset('window',34)
clf();
[solutionC]=CrankNicholson(temps,y0,G);
plot(temps,solutionC(1,:),'r')
plot(temps,solutionC(2,:),'g+')
plot(temps,solutionC(3,:),'b-')
xlabel('temps')
ylabel('évolution des trois populations')
xtitle('perturbation du cinquième état d''équilibre');

//courbes données par le solveur ODE

[solode]=ode(y0,0,temps,G);

xset('window',35)
clf();
plot(temps,solode(1,:),'r')
plot(temps,solode(2,:),'g+')
plot(temps,solode(3,:),'b-')
xlabel('temps')
ylabel('évolution des trois populations')
xtitle('perturbation du cinquième état d''équilibre donné par ode');
//courbes données par Euler modifié
[solutionEm]=solEuler1(temps,y0,G);
xset('window',36)
clf();

plot(temps,solutionEm(1,:),'r')
plot(temps,solutionEm(2,:),'g+')
plot(temps,solutionEm(3,:),'b-')
xlabel('temps')
ylabel('évolution des trois populations')
xtitle('perturbation du cinquième état d''équilibre donné par le schéma d''Euler modifié');

//courbes données par RK4
[solution]=solRK4(temps,y0,G);
xset('window',37)
clf();

plot(temps,solution(1,:),'r')
plot(temps,solution(2,:),'g+')
plot(temps,solution(3,:),'b-')
xlabel('temps')
ylabel('évolution des trois populations')
xtitle('perturbation du cinquième état d''équilibre donné par le schéma RK4');
