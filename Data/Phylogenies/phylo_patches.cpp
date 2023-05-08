// compila-se da seguinte forma c++ -O3 neutral_patch.cpp -o neutral_patch -lm -lgsl -lgslcblas
#include <math.h>
#include <sstream>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <fstream>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_double.h>
#include <unistd.h>
#include <algorithm>

using namespace std;
long idum;
gsl_rng *gerador;

struct grupo{
  int *size;
  int init_size;
  int K_init;
  double **mig;
  int ***ind_gen;
  int ***ind_gen_mig;
  int **filo;

};

struct rede{
  int Ngroups;
  int *C;
  int **Node;
};

struct populacao{
  double mut;
  double g;
  int ntraits;
  double migration;
  int diversity;
  int genoma;
  int name_max;
  int div_conf;

};

void initialization(struct populacao *popul, struct rede *network, struct grupo *group);
void migration(struct populacao *popul, struct rede *network, struct grupo *group, double mig);
void recombination(struct populacao *popul, struct rede *network, struct grupo *group);
void diversity(struct populacao *popul, struct rede *network, struct grupo *group, int t, int tmax, int conf, double mig, FILE *Temp, FILE *MS);
void definenet(struct rede *network);



int main(int ac, char **av)
{

  FILE *Temp, *MS;
  grupo group;
  rede network;
  populacao popul;

  int i, j, k, l, ano_mar, threshold, Isol, conf, config, t, tmax, cont, Pext, dham;
  double alt_mar, mig;

  long seed;  //define a seed

  char arq1[1000], arq2[1000];

  if (ac!=4)
    {
      cout  <<  "start the program like this:\n" << av[0]
            << "  [number_of_trees] [mean_migration_rate] [isolation_time]  \n"
            << endl;
      exit (-1);
    }

  /**** Leitura de dados por meio de um script ****/

  j = 0;
  config = atoi (av[++j]);
  popul.migration =  atof (av[++j]);
  Isol = atoi (av[++j]);


  cout  << "#invocation: ";
  for (int i=0; i<ac; i++){
    cout << av[i] << " ";
  }
  cout << endl;

  /// Parameters of the study
  network.Ngroups =  2;
  group.K_init = 200;
  group.init_size = 200;
  popul.genoma = 2000;
  popul.g = 0.9;
  popul.mut =  0.001;
  tmax = 2001;
  
  seed = time (NULL) * getpid();  //set the seed to system time
  gerador=gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(gerador, seed); //give the seed to random generator


  /// Time spent in connection for h = {0m, -10m, ..., -90m, -100m}
  double ts[11] = {1940,1820,1590,1430,1230,1010,820,620,480,290,140};

  double mig_ref = popul.migration;
  popul.migration = popul.migration*(ts[0]/ts[Isol/10]);


  network.C = new int[network.Ngroups];
  network.Node = new int*[network.Ngroups];
  for( i=0; i<network.Ngroups; i++ )
    network.Node[i] = new int[network.Ngroups];

  group.size = new int[network.Ngroups];                        /// Patches config


  group.ind_gen = new int**[network.Ngroups];
  for( i=0; i<network.Ngroups; i++ )
    group.ind_gen[i] = new int*[(int)(10*group.K_init)];
  for( i=0; i<network.Ngroups; i++ )
    for( j=0; j<(10*group.K_init); j++ )
      group.ind_gen[i][j] = new int[popul.genoma+1];

  group.ind_gen_mig = new int**[network.Ngroups];
  for( i=0; i<network.Ngroups; i++ )
    group.ind_gen_mig[i] = new int*[(int)(10*group.K_init)];
  for( i=0; i<network.Ngroups; i++ )
    for( j=0; j<(10*group.K_init); j++ )
      group.ind_gen_mig[i][j] = new int[popul.genoma+2];

  group.filo = new int*[1000];
    for( i=0; i<1000; i++ )
      group.filo[i] = new int[4];



    /// INPUT sea_level_data
    std::ifstream file("entrada_dados_nivel_mar.txt");
    std::vector<std::vector<double>> mar;
    std::string line;
    while (std::getline(file, line))
    {
      std::istringstream ss(line);
      std::vector<double> new_vec;
      double v;
      while (ss >> v)                 // populate the new vector
      {
        new_vec.push_back(v);
      }
      mar.push_back(new_vec);     // append it to the list of vectors
    }
    file.close();


    popul.div_conf=1;


  while( (++conf)<=config )
  {
    
      sprintf(arq1,"Pop%i_Isol%i_Mig%g_Temporal%i_.txt",group.K_init, Isol, mig_ref, conf);
      sprintf(arq2,"Pop%i_Isol%i_Mig%g_MS%i_.txt",group.K_init, Isol, mig_ref, conf);

      Temp = fopen(arq1,"w+");
      MS = fopen(arq2,"w+");
      fprintf(Temp,"pop rep mig tem esp abund1 abund2 \n");


      /// Initial conditions
      t = -1; ano_mar=0;
      threshold = -Isol;
      definenet(&network);
      initialization(&popul,&network,&group);


      while( (++t<tmax) )
  	  {
    
    	    recombination(&popul,&network,&group);
    
          if(t%10==0)
          {
            diversity(&popul,&network,&group,t,tmax,conf,mig,Temp,MS);
            ///cout<<t<<'\n';
          }
    
          if(t==mar[ano_mar][0])
          {
                alt_mar = mar[ano_mar][1];
                if(alt_mar>=threshold)
                {
                    mig = 0;
                }
                else
                {
                    mig=popul.migration;
                }
                ano_mar++;
          }
    
          if(Isol==0){mig=popul.migration;} /// Continuos migration (h>19.2m)
      	  migration(&popul,&network,&group,mig);
      
  	  }


      fclose(Temp);
      fclose(MS);
    }


}

void definenet( rede *network )
{
  int i, j;

  for( i=0; i<network->Ngroups; i++ )
    network->C[i] = 0;

  for( i=0; i<network->Ngroups; i++ )
    {
      for( j=0; j<network->Ngroups; j++ )
	if( i!=j )
	  {
	    network->Node[i][network->C[i]] = j;
	    network->C[i]++;
	  }
    }

}

void initialization(populacao *popul, rede *network, grupo *group)
{
  int i, j, k;
  for( i=0; i<network->Ngroups; i++ )
    {
      group->size[i] = group->init_size;


      for( j=0; j<group->size[i]; j++ )
        {
             for( k=0; k<popul->genoma; k++ )
		    {
                group->ind_gen[i][j][k] = 1;
		    }
		    group->ind_gen[i][j][popul->genoma] = 0;        /// Inicializando o nome x da espécie [x][]
        }

    }
}

void recombination(populacao *popul, rede *network, grupo *group)
{

  int s, i, j, k, partner, position, busca, aux_group_size, *aux_ind_gen, **aux_group_ind_gen;
  double recomb, dham, r;

  aux_ind_gen = new int[popul->genoma+1];

  aux_group_ind_gen = new int*[10*group->K_init];
  for( i=0; i<(10*group->K_init); i++ )
    aux_group_ind_gen[i] = new int[popul->genoma+1];

  for( i=0; i<network->Ngroups; i++ )
  {
          //  cout << group->size[i] << endl;
          aux_group_size = 0;
          if( group->size[i]==1 ){group->size[i] = 0;}


          for( s=0; s<group->K_init; s++ )
          {
  
            focal:
            j = (gsl_ran_flat(gerador, 0., 1.)*group->size[i]); /// random focal parent
            busca=0;
            similaridade:
            partner = (int)(gsl_ran_flat(gerador, 0., 1.)*group->size[i]);
            while( partner==j ) { partner = (int)(gsl_ran_flat(gerador, 0., 1.)*group->size[i]); }
            dham=0;
            for( k=0; k<popul->genoma; k++ ) { dham += group->ind_gen[i][j][k]*group->ind_gen[i][partner][k]; }
            if((dham/popul->genoma)<popul->g ){busca++; if(busca<=group->K_init) {goto similaridade;}else{goto focal;}} /// substituir K_init por groupsize
            if((dham/popul->genoma)>= popul->g )
            {
                for( k=0; k<popul->genoma; k++ )
                {
                    position = (int)(gsl_ran_flat(gerador, 0., 1.)*10);
                    if(position>=5){ aux_ind_gen[k] = group->ind_gen[i][j][k];} else{ aux_ind_gen[k] = group->ind_gen[i][partner][k]; ;}
  
                    r = gsl_ran_flat(gerador, 0., 1.); ///cout<<r<<'\n';
                    if( r<popul->mut )
                    {
                        if(aux_ind_gen [k]==1) { aux_ind_gen [k]=-1;} else aux_ind_gen [k]=1;   /// [+,-1]
                        ///if(aux_ind_gen [position]==1) { aux_ind_gen [position]--;} else aux_ind_gen [position]++;   /// [0,1]
                    }
  
                }
                aux_ind_gen [popul->genoma] = group->ind_gen[i][j][popul->genoma];
                for( k=0; k<popul->genoma+1; k++ ) { aux_group_ind_gen[aux_group_size][k] = aux_ind_gen[k]; }
                aux_group_size++;
            }
          }

        group->size[i] = aux_group_size;
        for( j=0; j<group->size[i]; j++ )
        {
          for( k=0; k<popul->genoma+1; k++ ) { group->ind_gen[i][j][k] = aux_group_ind_gen[j][k]; }
        }

    }

  for( i=0; i<(10*group->K_init); i++ )
    delete[] aux_group_ind_gen[i];
  delete[] aux_group_ind_gen;

  delete[] aux_ind_gen;

}

void migration(populacao *popul, rede *network, grupo *group, double mig)
{
  int NDinc[5001], NDdec[5001], k, i, j, ind_group;
  int num_gen[5001];
  double r;

  for( i=0; i<network->Ngroups; i++ )
    NDinc[i] = NDdec[i] = 0;

  for(i=0; i<5000; i++){num_gen[i]=0;}
  for( i=0; i<network->Ngroups; i++ )
  {

        for( j=0; j<group->size[i]; j++ )
        {
          r = gsl_ran_flat(gerador, 0., 1.);
          if( r<mig)                                          
            {

              ind_group = (int)(gsl_ran_flat(gerador, 0., 1.)*network->Ngroups);
              while(ind_group==i){ind_group = (int)(gsl_ran_flat(gerador, 0., 1.)*network->Ngroups); ;}

              for( k=0; k<popul->genoma+1; k++ )
              {
                  group->ind_gen_mig[ind_group][num_gen[ind_group]][k] = group->ind_gen[i][j][k];
              }
                num_gen[ind_group]++;

            }
            else
            {
              ind_group = i;
              for( k=0; k<popul->genoma+1; k++ )
              {
                  group->ind_gen_mig[ind_group][num_gen[ind_group]][k] = group->ind_gen[i][j][k];
              }
                num_gen[ind_group]++;
            }
        }
    }


    for( i=0; i<network->Ngroups; i++ )
    {
        group->size[i] = num_gen[i];
        for( j=0; j<group->size[i]; j++ )
        {
          for( k=0; k<popul->genoma+1; k++ ) { group->ind_gen[i][j][k] = group->ind_gen_mig[i][j][k]; }
        }

    }


}

void diversity(populacao *popul, rede *network, grupo *group, int t, int tmax, int conf, double mig, FILE *Temp, FILE *MS)
{

 int i, j, j_aux, k, k_aux, s, h, b, c, num_spec, num_gen, new_gen, equal_gen, cont, n, freq_total, name_max, aux_name_max, abundant, equa;
 double dham;
 int **gen, **freq, **sum, *freq_spec, *S, **spec, *contabilizado, *comu, **name;


 S = new int[(int)(10*group->K_init)];
 contabilizado = new int[(int)(10*group->K_init)];
 freq_spec = new int[(int)(10*group->K_init)];
 comu = new int[(int)(10*group->K_init)];


 for(i=0; i<(int)(10*group->K_init); i++) { contabilizado[i]=0; }
 for(i=0; i<(int)(10*group->K_init); i++) { freq_spec[i]=0; }
 for(i=0; i<(int)(10*group->K_init); i++) { comu[i]=1; }

   freq = new int*[10*group->K_init];
  for( i=0; i<(10*group->K_init); i++ )
    freq[i] = new int[network->Ngroups+1];

  for(i=0; i<(int)(10*group->K_init); i++)
  {
      freq[i][network->Ngroups]=1;
      for(j=0; j<network->Ngroups; j++){freq[i][j]=0; }
  } 

   gen = new int*[10*group->K_init];
  for( i=0; i<(10*group->K_init); i++ )
    gen[i] = new int[popul->genoma+1];

   spec = new int*[10*group->K_init];
  for( i=0; i<(10*group->K_init); i++ )
    spec[i] = new int[10*group->K_init];

   name = new int*[1000];
  for( i=0; i<1000; i++ )
    name[i] = new int[network->Ngroups+4];

  sum = new int*[1000];
  for( i=0; i<1000; i++ )
  {
    sum[i] = new int[network->Ngroups+1];
  }


  if(popul->div_conf==conf)
  {
    for( i=0; i<(10*group->K_init); i++ )
    {
        for( j=0; j<(10*group->K_init); j++ )
        {
            spec[i][j]=0;
        }
    }

    for( i=0; i<(10*group->K_init); i++ )
    {
        for( j=0; j<popul->genoma; j++ )
        {
            gen[i][j]=0;
        }
    }

    for( i=0; i<1000; i++ )
    {
        for( j=0; j<network->Ngroups+4; j++ )
        {
            name[i][j]=0;
        }
    }


    for( i=0; i<1000; i++ )
    {
        for( j=0; j<network->Ngroups+1; j++ )
        {
            sum[i][j]=0;
        }
    }

    for( i=0; i<(10*group->K_init); i++ )
    {
        S[i]=0;
    }

    group->filo[0][0] = 0;
    group->filo[0][1] = 0;
    group->filo[0][2] = 0;
    group->filo[0][3] = tmax;

    popul->name_max=1;
    popul->div_conf++;
  }


/// /// ///////////////////////////////////////////////////////////////////// GENOTYPE FREQUENCY [matrix] for large L

    /// Post Entries
    num_gen=0;
    for( i=0; i<network->Ngroups; i++ )
    {
        for( j=0; j<group->size[i]; j++ )
        {
            new_gen=0;
            for( k_aux=0; k_aux<num_gen; k_aux++ )
            {
                equal_gen=0;
                for( k=0; k<popul->genoma; k++ )
                {
                    if(gen[k_aux][k]==group->ind_gen[i][j][k]){ equal_gen++; } 

                }

                if(equal_gen==popul->genoma) 
                {
                    freq[k_aux][network->Ngroups]++; 
                    freq[k_aux][i]++;
                    goto same_gen;          
                }   else{new_gen++;}

            }

            if(new_gen==num_gen) 
            {
                for( k=0; k<popul->genoma+1; k++ )
                {
                    gen[num_gen][k] = group->ind_gen[i][j][k];
                }
                num_gen++;
                freq[k_aux][i]++;
            }
        same_gen:;
        }

    }



/// /// /////////////////////////////////////////////////////////////////////////////////// SPECIES RICHNESS

    spec[0][0]  = 0;
    contabilizado[0] = -1;
    cont = 1;
    comu [0] = 1;
    num_spec = 0;

    for(j=0; j<num_gen; j++)
    {



                for(i=0; i<cont; i++)              
                {
                    if(j==contabilizado[i])
                    {
                        goto j_next;
                    }
                }
                spec[num_spec][0] = j;
                contabilizado [0] = j;
        for( c=0; c<comu[num_spec]; c++ )
        {

            for( j_aux=0; j_aux<num_gen; j_aux++ )
            {
                h=0;
                for(i=0; i<cont; i++)             
                {
                    if(j_aux==contabilizado[i])
                    {
                        goto j_aux_next;
                    }
                h++;
                }




                if(h==cont)
                {
                    dham=0;
                    for( k=0; k<popul->genoma; k++ )
                    {
                        dham += gen[spec[num_spec][c]][k]*gen[j_aux][k];
                        
                    }
                        
                    if((dham/popul->genoma)>= popul->g )    
                    {
                        spec[num_spec][comu[num_spec]] = j_aux;
                        contabilizado [cont] = j_aux;

                        comu[num_spec]++;
                        cont++;

                    }
                }


            j_aux_next:;
            }

        }
                    num_spec++;

  j_next:;}



        freq_total = 0;
        for( j=0; j<num_spec; j++ )
        {
            for( k=0; k<comu[j]; k++ )
            {
                freq_spec[j] += freq[spec[j][k]][network->Ngroups];
            }
            freq_total+=freq_spec[j];
        }

        for( j=0; j<num_spec; j++ )
        {
            for( k=0; k<comu[j]; k++ )
            {
               S[spec[j][k]] = j;
            }

        }


/// /// /////////////////////////////////////////////////////////////////////////////////// NAMING SPECIES

        s=0;
        cont=0;
        aux_name_max = popul->name_max;
        while(cont<aux_name_max)
        {
            for( j=0; j<1000; j++ )
            {
                for(k=0; k<network->Ngroups+1; k++)
                {
                    sum[j][k] = 0;
                }

            }
            for( j=0; j<num_gen; j++ )
            {
                if(gen[j][popul->genoma]==cont)
                {
                    for(k=0; k<network->Ngroups+1; k++)
                    {
                        sum[S[j]][k] += freq[j][k];
                    }

                }
            }

            abundant=0;
            for( j=0; j<1000; j++ )             
            {
                if(sum[j][network->Ngroups]!=0)
                {
                        if(sum[j][network->Ngroups]>=abundant){ abundant = sum[j][network->Ngroups];}
                }
            }

    equa =0 ;
                                             
            for( j=0; j<1000; j++ )             
            {
                if(sum[j][network->Ngroups]!=0)
                {
                        name[s][0] = cont;                          
                        name[s][1] = j;                             
                        
                        for(k=3; k<network->Ngroups+4; k++)
                        {
                            name[s][k] = sum[j][k-3];
                        }

                        for(i=0; i<s; i++)                          
                        {
                            if((name[i][1]==name[s][1]) && (name[i][0]!=name[s][0]))
                            {
                                name[s][2] = name[i][0];
                                goto fusao;
                            }
                        }


                        if((name[s][5]==abundant) && (equa==0))
                        {
                            name[s][2] = cont ; equa=1;
                        } else {
                                    name[s][2] = popul->name_max;
                                    group->filo[popul->name_max][0] = name[s][2];
                                    group->filo[popul->name_max][1] = name[s][0];
                                    group->filo[popul->name_max][2] = t;
                                    group->filo[popul->name_max][3] = tmax;
                                    popul->name_max++;
                                    
                                }
                        
                    fusao:
                    s++;
                    
                }

            }


            cont++;
        } /// while             


        for(j=0; j<popul->name_max; j++)
        {
            cont=0;
            for(k=0; k<s; k++)
            {
                if(group->filo[j][0]==name[k][2]){ cont=1;}
            }
            if(group->filo[j][3]==tmax)
            {

                if(cont==0){ group->filo[j][3]=t;}
            }
        }



        for( i=0; i<s; i++ )
        {
            fprintf(Temp,"%i %i %g %i ", 2*group->K_init,1,mig,t);
            for( j=2; j<network->Ngroups+3; j++ )
            {
                if(j==2){fprintf(Temp,"%i ",name[i][j]+1);}
                else{fprintf(Temp,"%i ",name[i][j]);}
            }
            fprintf(Temp,"\n");
        }


        if(t==2000)
        {
            for( i=0; i<popul->name_max; i++ )
            {
                for( j=0; j<4; j++ )
                {
                    if(j==0 || j==1){fprintf(MS,"%i ",group->filo[i][j]+1);}
                    else {fprintf(MS,"%i ",group->filo[i][j]);}
                }
                fprintf(MS,"\n");
            }

        }


/// /// /////////////////////////////////////////////////////////////////////////////////// PASSING THE NAME

        for(i=0; i<s; i++)
        {
            for( j=0; j<num_gen; j++ )
            {
                if(name[i][1]==S[j])
                {
                    gen[j][popul->genoma+1] = name[i][2];
                }
            }
        }



for(j_aux=0; j_aux<num_gen; j_aux++)
{
    if(gen[j_aux][popul->genoma+1]!= gen[j_aux][popul->genoma])
    {
        for( i=0; i<network->Ngroups; i++ )
        {
            for( j=0; j<group->size[i]; j++ )
            {
                if(group->ind_gen[i][j][popul->genoma]==gen[j_aux][popul->genoma]) 
                {
                    equal_gen=0;
                    for( k=0; k<popul->genoma; k++ )
                    {
                        if(gen[j_aux][k]==group->ind_gen[i][j][k]){ equal_gen++; } 

                    }
                    if(equal_gen==popul->genoma) {group->ind_gen[i][j][popul->genoma] =  gen[j_aux][popul->genoma+1];}
                }
            }
        }
    }
}






    popul->diversity = num_spec;





  for( i=0; i<(10*group->K_init); i++ ){   delete[] gen[i];}
  delete[] gen;

  for( i=0; i<(10*group->K_init); i++ ){   delete[] spec[i];}
  delete[] spec;

  for( i=0; i<(10*group->K_init); i++ ){   delete[] freq[i];}
  delete[] freq;

  for( i=0; i<1000; i++ ){   delete[] name[i];}
  delete[] name;

  for( i=0; i<1000; i++ ){   delete[] sum[i];}
  delete[] sum;

  delete[] freq_spec;
  delete[] comu;
  delete[] contabilizado;
  delete[] S;

} 


