
printf('\n                   *********** Soyez les Bienvenus************                         \n');
printf('\n       Projet-Optimisation-Groupe6-KOUAGOU N. Jean et BALOGOUN O.A. Ismaïla         \n\n');
                 printf('\n\n                        Début                  \n\n');
//1.Lecture du fichier profil.txt
M=fscanfMat("C:\Users\Belong''s Telecoms\Desktop\"+"profil.txt");
x=M(:,1);
g=M(:,2);
printf('\n 1.Le fichier profil.txt est lu. Il a fourni les données recueillies par les géomètres\n Pour le profil du terrain, voir figure n°0\n');
//Affichage de la fonction g
xset('window',0)
clf()
plot(x,g)
xlabel('absisses des points de mesures')
ylabel('les altitudes de ces points')
xtitle('Profil du terrain')
//Affichage de la fonction g1
xset('window',1)
clf 
x1=x(1:10);
g1=g(1:10);
plot(x1,g1)
xlabel('absisses des points de mesures en mètre')
ylabel('les altitudes de ces points en mètre')
//2.pente locale maximale sur les dix points
x1=x(1:10);
g1=g(1:10);
p0=(g1(2)-g1(1))/(x1(2)-x1(1));
p=p0;
for i=3:10
   p1=abs((g1(i)-g1(i-1))/(x1(i)-x1(i-1)));
   p=[p,p1];
   
end
pmax1=max(p)//pente locale maximale sur les dix points
//printf('2.a) La pente locale maximale sur les dix premiers poinrs est pmax1=%f\n\n',pmax1);

//pente locale maximale sur tout le terrain
p0=(g(2)-g(1))/(x(2)-x(1));
p=p0;
for i=3:185
   p1=abs((g(i)-g(i-1))/(x(i)-x(i-1)));
   p=[p,p1];
   end
   
p=max(p)//pente locale maximale sur tout le terrain
printf('\n2. La pente locale maximale sur tout le terrain est p=%f\n',p);

//3. Définition de la matrice A,du vecteur b et de la matrice C
//defition de la fonction f qui prend en paramètre n et L et renvoie les matrice A, b et C

   funcprot(0);
function [A,b,C]=f(n,L)
    h=L/(n-1);
    B=ones(1,n);
    B(1,1)=1/2;
    B(1,$)=1/2;
    A=diag(h*B);
    b=alpha*h*ones(2*(n-1),1);
    a=eye(n-1,n);
    for i=1:(n-1)
        for j=1:n
            if (j==i)
                a(i,j)=-1;
            elseif (j-i==1)
                a(i,j)=1;
            end
        end
        
    end
    c=eye(n-1,n);
    for i=1:(n-1)
        for j=1:n
            if (j-i==1)
                c(i,j)=-1
            end
        end
        
    end
    C=[a;c];
endfunction
printf('\n 3.La matrice A,le vecteur b et la matrice C sont définis, ils permettent de faire tourner ce programme\n');
//4. Définition de la matrice Am1, inverse de A
printf('\n 4.La matrice Am1, l''inverse de A est définie\n');
   // Am1 est définie dans la fonction ci-apres (algorithme d'Uzawa), elle est donnée par Am1=inv(A)
// 5. Programmation de  l’algorithme d’Uzawa  
function [U,Z,k]=mimise_uzawa(A,b,C,U0,Z0,rho,eps,kmax,x)
    Am1=inv(A);
    U1=-0.5*Am1*C'*Z0+G;
    Z1=max(Z0+rho*(C*U1-b),0);
    Z=Z1;
    u=U0;
    U=U1;
    k=0;
    while norm(U-u)>=eps & k<=kmax
        u=U;
        U=-0.5*Am1*C'*Z+G;
        Z=max(Z+rho*(C*U-b),0);
        plot(x,U,'r');
        k=k+1;
        end
endfunction 
// définition de la fonction calculant la solution d'Uzawa mais qui n'affiche pas les graphes
function [U,Z,k]=mimise_uzawa_sans_figure(A,b,C,U0,Z0,rho,eps,kmax,x)
    Am1=inv(A);
    U1=-0.5*Am1*C'*Z0+G;
    Z1=max(Z0+rho*(C*U1-b),0);
    Z=Z1;
    u=U0;
    U=U1;
    k=0;
    while norm(U-u)>=eps & k<=kmax
        u=U;
        U=-0.5*Am1*C'*Z+G;
        Z=max(Z+rho*(C*U-b),0);
        k=k+1;
        end
endfunction 
printf('\n 5.L''algorithme d''Uzawa est programmé, il permet aussi de faire tourner ce programme\n')

//définition de la fonction J(u*,alpha) donnant le coût en fonction de alpha
function y=J(U,alpha)
    y=(U-G)'*A*(U-G);
endfunction

//5. (suite) Tracé de la solution superposée au profil du terrain à chaque itération
       
// Programme principal sur les dix points
L1=132.53498;
alpha=0.1;
n=10;
rho=2;
eps=0.1;
kmax=10;
G=g1;
Z0=ones(1,2*(n-1))';
U0=zeros(x1);

// Affichage de la solution d'Uzawa sur les dix points superposée à la courbe de g1 
xset('window',1)

[A,b,C]=f(n,L1);
[U,Z,k]=mimise_uzawa(A,b,C,U0,Z0,rho,eps,kmax,x1);
xlabel('abscisses des points de mesures')
ylabel('les altitudes de ces points')
xtitle('Tracé de la solution d''Uzawa superposée au profil du terrain à chaque itération')
printf('\n                      Veuillez patienter. . .      \n\n');

// Programme principal sur les 185 points
L=2709.6042;
alpha=0.1;
n=185;
rho=2;
eps=0.01;
kmax=20000;
G=g;
Z0=ones(1,2*(n-1))';
U0=zeros(x);
[A,b,C]=f(n,L);
//Affichage de la solution d'Uzawa superposée au profil du terrain à chaque itération

xset('window',2)
plot(x,g)
xlabel('abscisses des points de mesures')
ylabel('les altitudes de ces points')
xtitle('Tracé de la solution d''Uzawa superposée au profil du terrain à chaque itération')
[U,Z,k]=mimise_uzawa(A,b,C,U0,Z0,rho,eps,kmax,x);

printf('\n5.Pour le tracé de la solution superposée au profil du terrain à chaque itération, aller à la figure n°2 \n');


// 6. Tracé final superposé au profil du terrain
xset('window',3)
plot2d(x,[g U],[1 2])
legend('Profil du terrain','Profil final donné par Uzawa',[1,2],'ur')
printf('\n6.1 Pour le tracé final superposé au profil du terrain, aller à la figure n°3\n')
printf('6.2 Le nombre d''itérations nécessaires pour obtenir le profil de contrainte alpha=%.2f avec eps=%.2f est %i\n',alpha,eps,k);

// 7. Variation de la pente 
     for  alpha=[0.09,0.1,0.2,0.5,0.7];
     [A,b,C]=f(n,L);
    [U,Z,k]=mimise_uzawa_sans_figure(A,b,C,U0,Z0,rho,eps,kmax,x);
    printf('\nLe nombre d''itérations nécessaires pour alpha = %f est %i\n',alpha,k);
   end
   printf('On remarque que le nombre d''itérations diminue quand alpha augmente et augmente quand alpha diminue.\n\n');
//8. Vériﬁons que si on choisit alpha > p le programme s’arrête au bout d’une itération
for alpha =[0.62,0.63,0.65,0.7,0.9,1]
[A,b,C]=f(n,L);
[U,Z,k]=mimise_uzawa_sans_figure(A,b,C,U0,Z0,rho,eps,kmax,x);
   printf('\n.Le nombre d''itérations nécessaires pour alpha = %f est %i\n',alpha,k);
  //les contraintes en U* pour alpha>p
  V=C*U-b;
  printf('Les contraintes d''inégalités en U* pour alpha=%f sont:\n',alpha);
  disp(V);
   printf('On remarque que ces contraintes sont toutes inactives\n');
 end
 printf('\n Ainsi, toutes les contraintes sont inactives en U* pour alpha > p et le nombre d''itérations nécessaires est 1\n');
//9 Tracé de la courbe J(u*,α) en fonction de alpha
   xset('window',4)
     clf()
     yy=0;
for  alpha=[0.09,0.1,0.2,0.3,0.5,0.6];
     [A,b,C]=f(n,L);
  [U,Z,k]=mimise_uzawa_sans_figure(A,b,C,U0,Z0,rho,eps,kmax,x);
   y=J(U,alpha);
   yy=[yy,y];
//on ignorera le premier élément de yy pour la représention du graphe de J(u,alpha)
end
t_alpha=[0.09,0.1,0.2,0.3,0.5,0.6];//un vecteur contenant les valeurs de alpha
y=yy(2:$);//on a ignoré le premier élément de yy
plot(t_alpha,y);
xlabel('alpha')
ylabel('J(u*,alpha)')
xtitle('tracé de la courbe J(u*,α) en fonction de alpha')
printf('\n 9. Pour le tracé de la courbe J(u*,α) en fonction de alpha, aller à la figure n°4\n')
//10. Sur la solution obtenue a t-on retiré plus de terrain qu’on en a remis ou l’inverse?
//soit diff_u_g la différence '"des aires"' entre U* et G 
I=ones(n,1); 
//On a:
diff_u_g=(U-G)'*(A*I)// on trouve diff_u_g positif donc on a plus remis de terrain qu'on en a retiré
printf('\n 10. On trouve la différence '"des aires"' positive donc on a plus remis de terrain qu''on en a retiré\n');
  //Fin de la programmation
//Pour les questions 11 et 12 voir rapport

