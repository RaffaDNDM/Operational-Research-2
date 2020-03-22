#include <cplex.h>

/*********************************************
int TSPopt(instance *inst);
*********************************************/
int xpos(int i, int j, instance *inst) { return inst->xstart + i * inst->nnodes + j; }
/*********************************************/


int xpos(int i, int j, instance *inst) {

    if(i==j)
        error

    if(i>j)
        return xpos(j,i, inst)

    int pos = i * inst->nnodes +j; //se fosse piena
    int pos = i * inst->nnodes +j -(i+1)*(i+2)/2;

    return pos;
}

cplex
{
  //open cplex model (Problema vuoto)
  int error;
  CPXENVptr env = CPXopenCPLEX(&error);
  CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

  //CPXcreateprob restituisce errore, nome ambiente

    CPXLIBAPI
    CPXENVptr CPXPUBLIC
        CPXopenCPLEX (int *status_p);

    CPXLIBAPI
    CPXLPptr CPXPUBLIC
        CPXcreateprob (CPXCENVptr env, int *status_p, char const *probname_str);

    //costruzione modello popolato
    build_model(&inst, env, lp);

  //CPLEX's parameter setting

  CPXmipopt(env, lp); //alla fine la soluzione sarà in una sua struttura dati, la memorizza e poi gliela dobbiamo richiedere

  //use the optimal solution found by CPLEX

  //free and close cplex model
  CPXfreeprob(env, &lp); //deallocazione di puntatore
  //free and close cplex model
  CPXclose(&env); //deallocazione env
}


build_model()
{
    double zero = 0.0;
    char binary = 'B';

    char** cname = (char **) calloc(1, sizeof(char));
    cname[0] = (char *) calloc(100, sizeof(char));


    //di default cplex minimizza
    //add binary var xij for i<j
  	for ( int i = 0; i < inst->nnodes; i++ )
  	{
  		for ( int j = i+1; j < inst->nnodes; j++ )
  		{
  			sprintf(cname[0], "x(%d,%d)", i+1,j+1); // memorizzare nome xij solo per me (non interessa a c_plex) (scalato per utente di 1)
  			double obj = dist(i,j,inst); //costo di arco = distanza valore del costo nella funzione obiettivo
  			double lb = 0.0;
  			double ub = 1.0;

        //tipo di variabile binaria (no frazionaria, la devo passare a cplex con 'B' come valore)

        //env,lp=structs    1=#colonne obj=vettore di costi
        CPXnewcols(env, lp, #cols, vettore di costi, vettore di lb, vettore di ub, puntatore a 'B', array di stringhe) )

        if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) )
          print_error(" wrong CPXnewcols on x var.s");

        //passo indici (i,j) e restituisce posizione di xij nel modello
        //CPXgetnumcols(env,lp) = numero colonne in modello
        if ( CPXgetnumcols(env,lp)-1 != xpos(i,j, inst) ) //controllo che xij sia ultima colonna aggiunta=num_colonne-1 (perchè appena messa)
          print_error(" wrong position for x var.s");
  		}
  	}

    //add degree constraints
    // x-constraints (out- and in-degree at nodes)
    for ( int h = 0; h < inst->nnodes; h++ )  //degree constraints
    {
      int lastrow = CPXgetnumrows(env,lp); // numero righe attuali= indice riga che devo generare
      double rhs = 2.0; // termine noto del vincolo
      char sense = 'E'; //equality vincolo 'L' minore o uguale 'G' maggiore o uguale
      sprintf(cname[0], "degree(%d)", h+1); //nome vincolo

      //nuove righe solo una riga
      //CPXnewrows(env, lp, numero righe, vettore di termini noti, vettore di tipo di vincoli, NULL, cname)

      if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [x1]");
      for ( int i = 0; i < inst->nnodes; i++ )
      {
          if(i==h)
            //xhh no modello
            continue;

          //CPXchgcoef(env, lp, riga attuale, posizione di var , nuovo valore)
          //xpos gestisce sia (3,5) che (5,3)
        if ( CPXchgcoef(env, lp, lastrow, xpos(i,h, inst), 1.0) ) print_error(" wrong CPXchgcoef [x1]");
      }
    }

    if(verbose>=100)
      CPXwriteprob(env, lp, "model.lp, NULL");

    free(cname[0]);
    free(cname);
}
