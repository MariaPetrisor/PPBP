#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mersenne.c"

typedef struct nod
{
    int nr;
    double t_start, s_end;
    struct nod *urm;

} Burst;

Burst *adaug_sfarsit(Burst *ultim, Burst *nou) //Functia de adaugare la sfarsit de lista a unui nod nou
{
    if (ultim!=NULL)
        ultim->urm=nou;
    return nou;
}

Burst *creeaza_nod(double info_start, double info_end, int nrr) //Functia care creeaza un nod nou
{
    Burst *lista;
    if ((lista=(Burst *)malloc(sizeof(Burst)))==NULL)
    {
        printf("Eroare: memorie insuficienta\n");
        return NULL;
    }
    else
    {
        lista->t_start=info_start;
        lista->s_end=info_end;
        lista->urm=NULL;
        lista->nr=nrr;
        return lista;
    }
}

void listez(Burst *radacina) //Functia de listare a listei
{
    Burst *lista;
    for (lista=radacina; lista!=NULL; lista=lista->urm)
        printf("%d. Start:%.2lf End:%.2lf  \n", lista->nr, lista->t_start,lista->s_end);
}

void eroare() //Afiseaza mesaj de eroare in cazul in care alocarea dinamica nu este efectuata corect
{
    printf("Eroare alocare dinamica a memoriei. \n");
    exit(1);
}


void diviziune (double Nn, double *tau, double h) //Functia ce divizeaza intervalul [0,T] cu pasul h<beta
{
    int i;
    for(i=1; i<=Nn; i++)
        tau[i]=(double)i*h;
}


int simulare( double lambda, int T, double **t, double tauNn )//Functia de simulare a procesului Poisson de rata lambda prin care
                                                             //se genereaza momentele de inceput t[i] ale sesiunilor de burst
{
    double u, x;
    int  k = 0;
    long int i;

    (*t) = (double*) calloc (1,sizeof(double)); // t[k] = t[0] = 0;

    if ( !(*t) )
        eroare();
    while ( (*t)[k] < T && (*t)[k]<tauNn)
    {
        u = genrand_real1();
        x = (-1/lambda) * log(u);
        k++;
        (*t) = (double*) realloc ((*t), (k+1)*sizeof(double));
        if ( !(*t) )
            eroare();
        (*t)[k] = (*t)[k-1] + x;
    }

    if ( (*t)[k] > T )
        k--;

    return k;
}

double Pareto (double alfa, double beta) //Prin intermediul metodei inversarii simulam distributia Pareto, returnand durata sesiunilor
{
    double u, x, p;
    u = genrand_real2();
    p=1/(1-alfa);
    if(u<1-1/alfa)
        x=u*alfa*beta/(alfa-1);
    else
        x=beta*pow(alfa*(1-u),p);
    return x;

}


int main()
{
    time_t secunde;
    secunde = time(NULL);
    init_genrand(secunde);

    int T=1000, Nn, mediaD=2;
    double lambda = 5, beta, alfa=1.9, h=0.4;

    int n, i, j, k, m, nrr;
    double d,s;
    double *t, *tau;
    Burst *lista;
    Burst *radacina=NULL;
    Burst *nou, *ultim=NULL;

    int nr_ses=0, *nrP, nrt;
    int r=20000, l_packet=256;

    Nn=(floor)(T/h);

    beta=(float)(mediaD*(alfa-1))/alfa; //Calculul lui Beta

    int *B; //B[i] stoca numarul de sesiuni active la momentul tau[i]
    B = (int *)calloc((Nn+1), sizeof(int));
    if (!B)
        eroare();


    tau = (double *)calloc((Nn+1), sizeof(double)); //Notam cu tau[i] punctele de diviziune, tau[i]=i*h, unde h este pasul
    if (!tau)
        eroare();

    diviziune(Nn,tau,h);

    B[0]=13; // Am apelat in MATLAB comanda n=poissrnd(/\), /\=lambda*mediaD
    n=13;
    nrr=0;
    for(i=0; i<n; i++)
    {
        d=Pareto(alfa,beta);
        if(d>T) d=T;
        nrr++;
        nou=creeaza_nod(0,d,nrr);
        if (radacina==NULL)  radacina=nou;
        ultim=adaug_sfarsit(ultim, nou);
        j=0;
        k=(int)(d/h) + (fmod(d,h)?1:0);

        for(m=j; m<=k; m++)
            B[m]++;
    }

    nrt=simulare(lambda,T,&t,tau[Nn]); //Apelam functia de simulare care modifica prin referinta vectorul t,
                                       //atribuindu-i momentele de inceput ale sesiunilor de burst si, de asemenea,
                                       //memoram numarul de sesiuni generate, nrt

    for(i=0; i<nrt; i++)
    {
        d=Pareto(alfa,beta);
        s=t[i]+d; //momentul de sfarsit al sesiunii

        if(s>T) s=T; //daca momentul de sfarsit depaseste sfarsitul intervalului, T, va fi trunchiat la valoarea T
        nrr++;
        nou=creeaza_nod(t[i],s, nrr);
        if (radacina==NULL)  radacina=nou;
        ultim=adaug_sfarsit(ultim, nou);

        j=(int)(t[i]/h) + (fmod(t[i],h)? 1:0);
        k=(int)(s/h) + (fmod(s,h)? 1:0);

        for(m=j; m<=k; m++) //actualizam numarul de sesiuni de burst active in momentele din punctele de diviziune incluse
                            //in intervalul de burst tau[j]...tau[k]
            B[m]++;


    }
    nr_ses=nrt+n;

    nrP= (int *)malloc((nr_ses)* sizeof(int));
    if (!nrP)
        eroare();

    double t_packet= l_packet*8.0/r;

    lista=radacina;

    for(i=0; i<nr_ses; i++)
    {
        nrP[i]=(int)(((lista->s_end)-(lista->t_start))/t_packet);
        lista=lista->urm;

    }

    FILE *f1, *f2, *f3;

    f1 = fopen( "sesactive.txt", "w" );
    f2 = fopen( "sesactive1.txt", "w" );
    f3 = fopen( "pachete.txt", "w" );
    if ( f1 == NULL )
    {
        printf( "Eroare citire din fisierul sursa.\n" );
        exit(1);
    }
    if ( f2 == NULL )
    {
        printf( "Eroare citire din fisierul sursa.\n" );
        exit(1);
    }
    if ( f3 == NULL )
    {
        printf( "Eroare citire din fisierul sursa.\n" );
        exit(1);
    }
    for(i=0; i<=Nn; i++)
        fprintf(f1,"%d \n", B[i] ); //Scriem in fisierul sesactive.txt

    for(i=0; i<Nn; i++)
        fprintf(f2, "%lf %d \n",tau[i], B[i] ); //Scriem in fisierul sesactive1.txt

    for(i=0; i<nr_ses; i++)
        fprintf(f3, "%d \n", nrP[i] ); //Scriem in fisierul pachete.txt

    fclose(f1);
    fclose(f2);
    fclose(f3);


    free(t);
    free(B);
    free(tau);

    return 0;
}
