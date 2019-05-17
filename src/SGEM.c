#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

/* #define MATHLIB_STANDALONE */
/* #include <Rmath.h> */

/* compilation command: */
/* compile: gcc -o SGEM SGEM.c -Wall -I/usr/include -lm -lgsl -lgslcblas */
/* gcc main.c -o GEM -Wall -I/usr/share/R/include -lm -lRmath */

#define sqr(x)   ((x)*(x))

const gsl_rng *gBaseRand;       /* global rand number generator */

/* function declaration */
double walltime(double *t0);
int offset( int x, int y, int z , int xSize, int ySize) { return ( z * xSize * ySize ) + ( y * xSize ) + x ; }
double laplace(FILE *fp, double *Z,double *logL,int K,int nind,int nloci);/* estimate the marginal likelihood using Laplace approximation */
double aic(FILE *fp,double *Z,double *logL,int K,int nind,int nloci);
double bic(FILE *fp,double *Z,double *logL,int K,int nind,int nloci);
double marg_lik_est(FILE *fp,double *Z,double *logL,int K,int nind); /* approximate the marginal likelihood using Chibs method */
int *genotype_s(FILE *fp, int *id, int nloci, int nind, int line, int ploidy); /* extract genotype data, Structure/BAPS format */
int *genotype_f(FILE *fp2, FILE *fp5, int *nloci, int *nind, int ploidy);/* extract genotype data, FASTA format */
void showint(int *a, int r, int c);
void show(double *a, int r, int c);
int allele_count(int *X,int nloci,int nind,int *allcnt,int ploidy,int *Xd); /* count no of alleles J_l at each loci l */
/* int y_kl_1(double *Z,int k,int l,int i);  *//* No of allels observed for locus l in ind assigned to population k excluding indiviudual i */
int x_ilj(int *Xi, int i, int l, int j, int ploidy, int nind); /* No of copies of allele j at locus l in ind i */
double lik_ki(int k, int i, int L); /* likelihood score for ind i assigned to population k */
double *init_alloc(const gsl_rng *gBaseRand,int K,int nind,int nindK,double c, int c0); /* initialize allocation vector prior to EM-algorithm */
int find_max(int *allcnt,int nloci);
void saveprint(FILE *fp,double *Z, int *id, int K, int nind);
int sum_alleles(int *allcnt,int nloci);
double update_contribution(double Z_ik,int *Xd,int i,int l,int j,int ploidy, int nind);
double *s_step(const gsl_rng *gBaseRand,double *Z,int nind,int K,int nindK,int t);
double y_klj(double *Z,int k,int nind, int *Xd, int l, int j, int ploidy);
void save_alloc(FILE *fp,double *Z,double *Zm,int nind,int K);
void recalc_mean(double *Zm,int nind,int K,int nit);

int main(int argc, char **argv){

  FILE *fp1,*fp2,*fp5,*fps,*fpZm,*fprho,*fpr,*fpm,*fp6;
  char *inp,data[100],zinit[200];
  int z0=0,ntot,nmax,ploidy,ncnt=0,K,line,nind,xilj,nloci,verbosity,k0=1,xil,conv,form,con;
  int git,i,j,k,l,m=1,nit,burnin,nmc,index,t0=0,i1,knindi,knindi1;//,flag_l;
  int *X,*allcnt,*Xd,*id;
  double *Z,*Zm,*rho,*rhom,*yklj,*yklj2,ykl,gamma_lj=1.,gamma_l0,li_ik,up_j,lo_j,pDat2,*logL,*L,pl,lpl,po,tol,diff,tmpn,la,a,b,tmp_y=0.,logLt,logLmax, logLt_old = 0.;
  float tol2,tmpf;
  double startTime, elapsedTime;
  double clockZero = 0.0;
  unsigned long randSeed;
  /* specifying to use Mersenne twister MT-19937 as the uniform PRNG */
  gBaseRand = gsl_rng_alloc(gsl_rng_mt19937);
  randSeed = rand();                    /* returns a non-negative integer */
  srand(randSeed);                    /* initialization for rand() */

  gsl_rng_set (gBaseRand, randSeed);    /* seed the PRNG */
  gsl_rng_env_setup();


  startTime = walltime(&clockZero);
  /******** read data from input file ********/

  if(argc > 1)
    inp = argv[1];
  else
    inp = "input.txt";

  if(!(fp1 = fopen(inp, "r"))){
    printf("Error, could not open user defined parameter file!\n");
    return 0;
  }
  fscanf(fp1," %*[^\n]"); /* skip header in file input.txt */
  fscanf(fp1,"%d%*[^\n] %d%*[^\n] %d%*[^\n] %d%*[^\n] %d%*[^\n] %s%*[^\n] %d%*[^\n] %d%*[^\n] %f%*[^\n] %d%*[^\n] %d%*[^\n] %d%*[^\n] %s%*[^\n] %d%*[^\n]",&ploidy,&form,&line,&nind,&nloci,data,&verbosity,&K,&tol2,&nit,&burnin,&git,zinit,&con);
  fclose(fp1);

  tol = (double)tol2;
  conv = nind*K;
  pl = (double)1./K; /* prior probability for inverse lambda, assuming a flat distribution */
  lpl = log(pl);
  nmc = nit + burnin;
  if(!(id = (int *)malloc(nind*sizeof(int)))){ 
    printf("Error in allocating memory in main program\n\n"); 
    return 0;
  }
  if(!(fp2 = fopen(data, "r"))){
    printf("Error, could not open data file!\n");
    return 0;
  }
 
  if(form==1) X = genotype_s(fp2, id, nloci, nind, line, ploidy);
  else{ 
    fp5 = fopen(data, "r");
    X = genotype_f(fp2, fp5, &nloci, &nind, ploidy);
    fclose(fp5);
    /* printf("**************\nNB! nloci = %d, nind = %d\n*************\n",nloci,nind); */
  }
  fclose(fp2);

  if(!(allcnt = (int *)malloc(nloci*sizeof(int))) ||
     !(Xd = (int *)malloc(nind*ploidy*nloci*sizeof(int))) ||
     !(L = (double *)malloc(K*sizeof(double))) ||
     !(logL = (double *)malloc(conv*sizeof(double)))){ 
    printf("Error in allocating memory in main program\n\n"); 
    return 0;
  }
  xil = ploidy;
  
  if(!allele_count(X,nloci,nind,allcnt,ploidy,Xd)){ /* count alleles per locus */
    /* Xr = reduce(X,nloci,nind,ploidy,allcnt);  *//* remove fixed loci */
    if(verbosity==1) printf("Need to remove loci with fixed alleles to improve computational efficiency!\n");
  }

  /* printf("Xd:\n"); */
  /* showint(Xd,nind,ploidy*nloci); */

  nmax = find_max(allcnt,nloci); /* find maximum number of alleles at any loci */
  ntot = sum_alleles(allcnt,nloci); /* sum total number of alleles at all loci */
  if(verbosity==1) printf("maximum number of alleles at any loci: %d\nTotal number of unique alleles in population: %d\n",nmax,ntot);
 
  if(git==0)
  // step 0: randomise starting allocation vector
    Z = init_alloc(gBaseRand, K, nind,conv,pl,1);
  else{
    if(!(Z = (double *)malloc(conv*sizeof(double)))){ 
      printf("Error in allocating memory in main program\n\n"); 
      return 0;
    }
    if(!(fp6 = fopen(zinit, "r"))){
      printf("Error, could not open initial allocation file: %s\n",zinit);
      return 0;
    }
    for(i = 0; i < nind; i++)
      for(j = 0; j < K; j++){
	fscanf(fp6,"%f",&tmpf);
	Z[j*nind+i] = (double)tmpf;
      }
    fclose(fp6);
  }
  /* printf("Z_0:\n"); */
  /* show(Z,nind,K); */
  
  if(!(yklj = (double *)malloc(nmax*sizeof(double))) ||
     !(rho = (double *)calloc(conv,sizeof(double))) ||
     !(rhom = (double *)calloc(conv,sizeof(double))) ||
     !(Zm = (double *)calloc(conv,sizeof(double))) ||
     /* !(Zs = (double *)calloc(conv,sizeof(double))) || */
     !(yklj2 = (double *)malloc(K*nmax*nloci*sizeof(double)))){
    printf("Error in allocating memory in main program\n\n");
    return 0;
  }
  if(!(fps = fopen("allocations.txt", "w")) ||
     !(fpr = fopen("mixture_weights.txt", "w"))){
    printf("Error, could not open data file!\n");
    return 0;
  }
  // step 1: start main loop GEM-iterations until convergence
  /* for(m = 0; m < nmc; m++){ */
  while(m <= nmc){
  /* while(it < nit){ */
  /* while(m < 5){ */
  /* m = 0; */
    logLt = 0.;
    if(burnin > 0) printf("GEM algorithm iteration number: %d\nE-step...\n",m);
    else printf("GEM algorithm iteration number: %d\nE-step...\n",k0);
    ncnt = 0;
    // step 1.1: loop over individuals 1...nind
    for(i = 0; i < nind; i++){ // i < nind
      if(i > 0) i1 = i - 1;
      else i1 = nind - 1;
      tmpn = -1.0E6;
      // step 1.1.1: loop over populations 1...K
      for(k = 0; k < K; k++){
	knindi = k*nind+i;
	knindi1 = k*nind+i1;
    	// loop over locus
    	li_ik = 0.;


    	for(l = 0; l < nloci; l++){
    	  // loop over alleles at loci l
    	  ykl = 0.;
    	  gamma_l0 = allcnt[l]*gamma_lj;
	  /* flag_l = 1; */
	  /* for(i1 = 0; i1 < ploidy; i1++) */
	  /*   if(Xd[(ploidy*l+i1)*nind+i] == -9) flag_l = 0; */
	  /* if(flag_l==1){ */
	    // for(j = 0; j < allcnt[l]; j++){
	    //yklj[j] = y_klj_1(Z,k,i,nind, Xd, l, j, ploidy);
	    //ykl += yklj[j];
	    //}
	  for(j = 0; j < allcnt[l]; j++){
	    //printf("hit1\n");
	    index = offset(k,l,j,K,nloci);
	    //printf("index = %d\n",index);
	    if(i == 0){
	      //printf("hit2\n");
	      yklj2[index] = y_klj(Z,k,nind, Xd, l, j, ploidy); // calculate contribution of all ind in pop k, locus l, allele j
	      /* printf("Prior removing contr of ind 1, pop %d, loci %d, allele %d: yklj2[%d] = %f\n",k+1,l+1,j+1,index,yklj2[index]); */
	    }
	    else tmp_y = update_contribution(Z[knindi1],Xd,i1,l,j,ploidy,nind); //-uppdate_contribution(Z_old[knindi1],Xd,i1,l,j,ploidy,nind);//offset(k,l,j,K,nloci) // change the contribution of individual i-1 (because of updated allocation Z[i-1,k])
	      // remove contribution from ind i
	    //printf("hit3\ni = %d, j = %d, k = %d, l = %d\nMemory allocated: %d\n",i,j,k,l,K*nmax*nloci);
	    /* if(i==1) printf("Before: yklj2[%d] = %f, tmp_y = %f, Z[%d] = %f, i-1 = %d\n",index,yklj2[index],tmp_y,knindi1,Z[knindi1],i1); */
	    yklj2[index] += tmp_y - update_contribution(Z[knindi],Xd,i,l,j,ploidy,nind); // remove contribution from ind i
	    //printf("hit3.5\n");
	    /* if(i==1) printf("After removing contr of ind 1: yklj2[%d] = %f\n",index,yklj2[index]); */
	    ykl += yklj2[index];
	    /* printf("loci = %d, allele nr = %d, yklj2[%d] = %f, tmp_y = %f, ",l+1,j+1,index,yklj2[index],tmp_y); */
	    //printf("hit4\n");
	      // step 1.1.1.1: calculate x and y conditional on current allocation z_ik
	    xilj = x_ilj(Xd, i, l, j, ploidy, nind);
	    
	      // step 1.1.1.2: calculate the current likelihood score
	      /* if(k==0 && l==0) printf("xilj: %d\ti: %d\tk: %d\tl: %d\tj: %d\tallcnt: %d\tykl: %f\tyklj[%d]: %f\n",xilj,i,k,l,j,allcnt[l],ykl,j,yklj[j]); */
	    lo_j = yklj2[index] + gamma_lj;
	    up_j = (double)xilj + lo_j;
	    li_ik += lgamma(up_j)-lgamma(lo_j); //gamma(up_j)/gamma(lo_j);
	    /* if(k==0 && l==0) printf("i: %d\tk: %d\tl: %d\tj: %d\tgamma_u[%f]: %f\tgamma_l[%f]: %f\tLikelihood: %f\n",i,k,l,j,up_j,lgamma(up_j),lo_j,lgamma(lo_j),li_ik); */
	      /* if(g_l==0.) printf("NB! gamma lower == 0\ti: %d\tk: %d\tl: %d\tj: %d\tyklj: %f\tgamma_lj: %f\tlo_j: %f\n",i,k,l,j,yklj[j],gamma_lj,lo_j); */
	    /* } */
	  }
	  up_j = ykl+gamma_l0;
	  lo_j = (double)xil + up_j;
	  li_ik += lgamma(up_j)-lgamma(lo_j);
	  /* if(i==1) printf("ind: %d, pop: %d, loci: %d, ykl: %f, xil: %d, up_j: %f, lo_j: %f, li_ik: %f\n",i+1,k+1,l+1,ykl,xil,up_j,lo_j,li_ik); */
	  /* printf("ind no: %d pop no: %d loci no %d: up_j: %f, lo_j: %f, gamma_u: %f, gamma_l: %f, lik: %f\n",i+1,k+1,l+1,up_j,lo_j,lgamma(up_j),lgamma(lo_j),li_ik); */
	  
	}
	  /* printf("\n"); */
	logL[knindi] = li_ik;
  	/* if (verbosity==1) printf("***** knindi = %d ******\nLogLikelihood for ind %d at pop no %d: %f\n",knindi,i+1,k+1,logL[knindi]); */
  	if(tmpn < logL[knindi] && logL[knindi] < 1000000.) /* find maximum logLikelihood */
  	  tmpn = logL[knindi];
	

    	/* L[k] = li_ik/K; */
    	/* LikPri += L[k]; */ // marginal likelihood
      }
      // step 1.1.1.3: calculate new posterior probability p(z_ik|x,y)
      
      /* pDat = 0.; */
      pDat2 = 0.;
      logLmax = -1.0E10;
      /* printf("tmpn = %f\n",tmpn); */
      for(k = 0; k < K; k++){
	knindi = k*nind+i;
  	if(logL[knindi] > -1000000. && logL[knindi] < 1000000.) L[k] = exp(logL[knindi]+lpl-tmpn);
  	else L[k] = 0.;
  	/* pDat += L[k]; */
  	pDat2 += L[k];
	if(logL[knindi] > logLmax) logLmax = logL[knindi];
	/* printf("L[%d] = %f, logL[%d] = %f, lpl = %f, tmpn = %f\n",k,L[k],knindi,logL[knindi],lpl,tmpn); */
      }
      if(verbosity==1) printf("Likelihood of ind %d belonging to alternative populations:\n",i+1);
      if(verbosity==1) show(L,1,K);
      if(verbosity==1) printf("Posterior probabilities for induvidual %d\n",i+1);
      for(k = 0; k < K; k++){
	knindi = k*nind+i;
  	po =  L[k]/pDat2;
	if(con == 1){
	  diff = sqr(po-rho[knindi]);
	  if(diff/po < tol) ncnt++;
	}
    	rho[knindi] = po;
    	if(verbosity==1) printf("%f ",po);
      }
      if(verbosity==1) printf("\n");
      logLt += logLmax;
    }
    printf("max loglikelihood: %f\n",logLt);
    if(con==1) printf("number of converged allocations: %d/%d\n",ncnt,conv);
    if(con==1 && burnin==0 && ncnt==conv) t0 = 1;
 
    /* if(burnin > 0 || t0 == 1) m++; */
   
    // step 1.2: calculate new convergence score after updating all posterior probabilities for z
    printf("E-step done\nS-step...\n");
    if(verbosity==1) printf("Draw indicator variable z_i for each individual i from the multinomial distribution z_i~multinom(1,rho)\n");
    /* printf("old allocation matrix:\n"); */
    /* show(Z,nind,k); */
    /* printf("mixture proportion matrix:\n"); */
    /* show(rho,nind,k); */

    Z = s_step(gBaseRand,rho,nind,K,conv,1);


    if(con == 2 && burnin==0){
      diff = sqr(logLt-logLt_old);
      printf("old logL: %f, new logL: %f, diff: %f, tol: %f\n",logLt_old,logLt,diff/logLt,tol);
      logLt_old = logLt;
      if(fabs(diff/logLt) < tol) t0 = 1;
      /* else t0 = 0; */
    }
    if(con == 3 && burnin==0){
      diff = fabs(logLt-logLt_old);
      printf("old logL: %f, new logL: %f, diff: %f, tol: %f\n",logLt_old,logLt,diff,tol);
      logLt_old = logLt;
      if(diff < tol) t0 = 1;
    }
    if((m > burnin && burnin > 0) || t0==1){ 
      save_alloc(fps,Z,Zm,nind,K);
      save_alloc(fpr,rho,rhom,nind,K);
      z0++;
      if(burnin == 0) m++;
      printf("save allocations, iteration no %d\n",z0);
    }
    if(burnin > 0) m++;
    else k0++;
    /* printf("new allocation matrix:\n"); */
    /* show(Z,nind,k); */
  }
  elapsedTime = walltime(&startTime)/60.0;
  free(yklj2);
  fclose(fps);
  fclose(fpr);
  printf("Elapsed time from start (in minutes): %f\n",elapsedTime);
  recalc_mean(Zm,nind,K,nit);
  recalc_mean(rhom,nind,K,nit);
  if(!(fpm = fopen("model_selection.txt", "w"))){
    printf("Error, could not open model selection save file!\n");
    return 0;
  }
  /* ma = marg_lik_est(fpm,Zm,logL,K,nind); */
  /* printf("marginal likelihood: %f\n",ma); */
  la = laplace(fpm,Z,logL,K,nind,nloci);
  printf("Laplace approximation: %f\n",la);
  a = aic(fpm,Zm,logL,K,nind,nloci);
  printf("AIC: %f\n",a);
  b = bic(fpm,Zm,logL,K,nind,nloci);
  fclose(fpm);
  printf("BIC: %f\n",b);
  if(!(fpZm = fopen("mean_allocation.txt", "w")) ||
     !(fprho = fopen("mean_mixture_proportion.txt", "w"))){
    printf("Error, could not save allocation files!\n");
    return 0;
  }
  /* save mean allocations and mixture proportions to files */   
  saveprint(fpZm,Zm,id,K,nind);
  fclose(fpZm);
  saveprint(fprho,rhom,id,K,nind);
  fclose(fprho);

  free(id);
  free(L);
  free(logL);
  free(allcnt);
  free(Xd);
  free(Z);
  free(yklj);
  free(rho);
  free(rhom);
  free(Zm);
  /* free(Zs); */

  return 1;

}

double laplace(FILE *fp, double *Z,double *logL,int K,int nind,int nloci){/* estimate the marginal likelihood using Laplace approximation */

  double a,pi_k,tmp;
  int i,j,jmax;
  double pmax,logLt=0.,detH=0.,pm=0.,pr;

  pr = (double) log(1./K);

  for(i = 0; i < nind; i++){
    pmax = 0.;
    for(j = 0; j < K; j++)
      if(Z[j*nind+i] > pmax){
	pmax = Z[j*nind+i];
	jmax = j;
      }
    logLt += Z[jmax*nind+i]*logL[jmax*nind+i]; //*log(Z[jmax*nind+i]); // likelihood
    pm += Z[jmax*nind+i]*pr; // prior
  }
  for(j = 0; j < K; j++){ 
    tmp = 0.;
    for(i = 0; i < nind; i++)
      tmp += Z[j*nind+i]; 
    pi_k = (double)tmp/nind;
    detH += log(tmp) - 2.*log(pi_k);
  } 

  a = (double) (logLt + pm + K*log(2.*M_PI)/2.-detH/2.);
  fprintf(fp,"************* Laplace approximation ***********\n************** following Holmes et al. (2012), PLoS One ********************\n");
  fprintf(fp,"log(p(data|model)): %f, loglikelihood: %f, K: %d, 2*K*log(2.*M_PI)/2.: %f, det(Hessian)/2: %f\n\n",a,logLt,K,(double)K*(nind+1.)*log(2.*M_PI)/2.,detH/2.);


  return a;

}
void recalc_mean(double *Zm,int nind,int K,int nit){
  
  int i,k;
 
  for(i = 0; i < nind; i++){
    for(k = 0; k < K; k++){
      Zm[k*nind+i] = (double)Zm[k*nind+i]/nit;
    }
  }

}
void save_alloc(FILE *fp,double *Z,double *Zm,int nind,int K){

  int i,k;
 
  for(i = 0; i < nind; i++){
    for(k = 0; k < K; k++){
      if(k<(K-1))
	fprintf(fp,"%f ",Z[k*nind+i]); // save allocations to file
      else
	fprintf(fp,"%f",Z[k*nind+i]);
      Zm[k*nind+i] += Z[k*nind+i];
    }
    fprintf(fp,"\n"); 
  }

}
double *s_step(const gsl_rng *gBaseRand,double *Z,int nind,int K,int nindK,int t){
  // generate new allocation matrix Z = [z_1,z_2,....,z_n]
  int i,k;
  unsigned int *tmp;
  double *Zn,*pr;

  if(!(tmp = (unsigned int *) malloc(K*sizeof(int))) ||
     !(pr = (double *) malloc(K*sizeof(double))) ||
     !(Zn = (double *) malloc(nindK*sizeof(double)))){ 
    printf("Error in allocating memory in function s_step\n\n"); 
    return NULL;
  }

  for(i = 0; i < nind; i++){
    for(k = 0; k < K; k++){ 
      pr[k] = Z[k*nind+i]; // copy rho_k, k=1,...,K to temp vector
      tmp[k] = 0; // reset output tmp vector
    }
    /* sample from multinomial distibution to generate allocations for all individuals */
    /* gsl_ran_multinomial (gBaseRand, nGenotype, popsize_i, pG, numG); */
    gsl_ran_multinomial (gBaseRand, K, t, pr, tmp);
    /* rmultinom(t, pr, K, tmp); */
    for(k = 0; k < K; k++) Zn[k*nind+i] = (double)tmp[k]/t; // copy result to updated indicator matrix
  }

  free(tmp);
  free(pr);

  return Zn;

}
double update_contribution(double Z_ik,int *Xd,int i,int l,int j,int ploidy, int nind){

  double cn;
  int xi1;

  xi1 = x_ilj(Xd, i, l, j, ploidy, nind);
  cn = (double)Z_ik*xi1;

  return cn;
}
double walltime( double *t0 )
{

  double mic, time;
  double mega = 0.000001;
  struct timeval tp;
  struct timezone tzp;
  static long base_sec = 0;
  static long base_usec = 0;

  (void) gettimeofday(&tp,&tzp);
  if (base_sec == 0)
    {
      base_sec = tp.tv_sec;
      base_usec = tp.tv_usec;
    }

  time = (double) (tp.tv_sec - base_sec);
  mic = (double) (tp.tv_usec - base_usec);
  time = (time + mic * mega) - *t0;
  return(time);
}
double y_klj(double *Z,int k,int nind, int *Xd, int l, int j, int ploidy){

 /* No of copies of allele j at locus l in ind assigned to pop k */

  double cnt=0.;
  int xi1,i1;

  for(i1 = 0; i1 < nind; i1++){
    xi1 = x_ilj(Xd, i1, l, j, ploidy, nind);
    cnt += (double)Z[k*nind+i1]*xi1;
  }

  return cnt;

}
int sum_alleles(int *allcnt,int nloci){

  int ntot = 0,i;

  for(i = 0; i < nloci; i++) ntot += allcnt[i];

  return ntot;

}
double bic(FILE *fp,double *Z,double *logL,int K,int nind,int nloci){

  double a;
  int i,j,jmax;
  double pmax,logLt=0.;


  for(i = 0; i < nind; i++){
    pmax = 0.;
    for(j = 0; j < K; j++)
      if(Z[j*nind+i] > pmax){
	pmax = Z[j*nind+i];
	jmax = j;
      }
    logLt += Z[jmax*nind+i]*logL[jmax*nind+i]; // likelihood
 
  }

  a = -2.*logLt + nind*K*log(nind);
  fprintf(fp,"************* BIC ***********\n");
  fprintf(fp,"BIC: %f, loglikelihood: %f, number of parameters: %d, K: %d\n\n",a,logLt,nind*K,K);

  return a;

}
double aic(FILE *fp,double *Z,double *logL,int K,int nind,int nloci){

  double a;
  int i,j,jmax;
  double pmax,logLt=0.;


  for(i = 0; i < nind; i++){
    pmax = 0.;
    for(j = 0; j < K; j++)
      if(Z[j*nind+i] > pmax){
	pmax = Z[j*nind+i];
	jmax = j;
      }
    logLt += Z[jmax*nind+i]*logL[jmax*nind+i]; // likelihood
 
  }

  a = -2.*logLt + 2.*nind*K;
  fprintf(fp,"************* AIC ***********\n");
  fprintf(fp,"AIC: %f, loglikelihood: %f, number of parameters: %d, K: %d\n\n",a,logLt,nind*K,K);
 


  return a;

}
double marg_lik_est(FILE *fp,double *Z,double *logL,int K,int nind){

  int i,j,jmax;
  double ma=0.,pmax,logLt=0.,zm=0.,pm=0.,pr;

  pr = (double) log(1./K);

  for(i = 0; i < nind; i++){
    pmax = 0.;
    for(j = 0; j < K; j++)
      if(Z[j*nind+i] > pmax){
	pmax = Z[j*nind+i];
	jmax = j;
      }
    logLt += Z[jmax*nind+i]*logL[jmax*nind+i]; // likelihood
    zm += Z[jmax*nind+i]*log(Z[jmax*nind+i]); // posterior
    pm += Z[jmax*nind+i]*pr; // prior
  }

  ma = logLt + pm - zm;
  fprintf(fp,"************* Marginal likelihood estimation ***********\n");
  fprintf(fp,"p(data|model): %f, loglikelihood: %f, logprior: %f, logposterior: %f, K: %d\n\n",ma,logLt,pm,zm,K);
 
  return ma;

}
void saveprint(FILE *fp,double *Z, int *id, int K, int nind){

  //FILE *fp;
  int i,j;

  // if(!(fp = fopen("mean_alloc.out", "w"))){
  // printf("Error, could not open save file in function saveprint\n");
  //  exit(-1);
  //}
  fprintf(fp,"  ");
  for(j = 0; j < K; j++)
    fprintf(fp,"%d ",j+1);
  fprintf(fp,"\n");

  for(i = 0; i < nind; i++){
    fprintf(fp,"%d ",id[i]);
    for(j = 0; j < K; j++)
      fprintf(fp,"%f ",Z[j*nind+i]);
    fprintf(fp,"\n");
  }
  
  //fclose(fp);

}
int *genotype_f(FILE *fp2, FILE *fp5, int *nloci, int *nind, int ploidy){

  int i,j,n2,*L,tmpnind,tmploci,n;
  char c,c1;
 
  i = 0;
  j = 0;
  while (c != EOF){
    c = fgetc(fp5);
    /* if(c == '>') printf("new induvidual %d...\n",i+1); */
    /* if(i%2 != 0){ */
 
    /*  /\* NB! MISSING DATA?? *\/ */
    /* } */
    if(c == '\n'){ 
      i++; // count rows in data file
      if(i%2 == 0){ 
	*nloci = j;
	tmploci = j;
      }
      j = 0; 
      //printf("hit\n");
    }
    else j++; // counts column in data file (base positions)
  }

  *nind = floor(i/2);
  tmpnind = floor(i/2);

  if(ploidy==2){ 
    *nloci /= 2;
    tmploci /= 2;
  }

  n2 = ploidy*tmploci;

  if(!(L = (int *) calloc(tmpnind*n2,sizeof(int)))){ 
    printf("Error in allocating memory in function genotype_f\n\n"); 
    return NULL;
  }

  i = 0;
  j = 0;
  n = 0;
  while (c1 != EOF){
    c1 = fgetc(fp2);
    if(i%2 != 0){
      /* if(c1 == 'A' || c1 == 'a'){  */
      /* 	L[j*tmpnind+i] = 0; */
      /* 	printf("hit1: ind %d have an 'A' at position %d\n",i/2,j); */

      /* } */
      if(c1 == 'T' || c1 == 't') L[j*tmpnind+n] = 1;
	/* printf("ind %d, pos %d: 'T'\n",i/2+1,n+1); */
      
      else if(c1 == 'G' || c1 == 'g') L[j*tmpnind+n] = 2;
      else if(c1 == 'C' || c1 == 'c') L[j*tmpnind+n] = 3;
      else if(c1 == 'N' || c1 == 'n' || c1 == '-') L[j*tmpnind+n] = -9; /* missing data */
   
    }
    if(c1 == '\n'){ 
      i++; // count rows in data file
      j = 0; 
      if(i%2 != 0) n = floor(i/2);
	/* printf("%d ",n); */
      /* printf("hit: new line %d\n",i); */
    }
    else j++; // counts column in data file (base position)
  }

  /* showint(L,tmpnind,tmploci); */
  /* printf("tmpnind = %d, tmploci = %d\n",tmpnind,tmploci); */
  return L;

}
int *genotype_s(FILE *fp, int *id, int nloci, int nind, int line, int ploidy){

  int k,i,j,n2 = 2*nloci, ni2 = 2*nind,*L;

  if(ploidy==2){
    /* allocating memory for variables used in function genotype_s */

    if(!(L = (int *) malloc(nind*n2*sizeof(int)))){ 
      printf("Error in allocating memory in function genotype_s\n\n"); 
      return NULL;
    }
 
    if(line == 1){
      for(i = 0; i < nind; i++){
	fscanf(fp,"%d",&id[i]); // id nr
	for(j = 0; j < n2; j++)
	  fscanf(fp,"%d",&L[j*nind+i]); 
	fscanf(fp,"\n");
      }
    }
    else{
      for(i = 0; i < ni2; i++){
	k = floor(i/2);
	fscanf(fp,"%d",&id[k]); // id nr
	for(j = 0; j < nloci; j++){
	  if(i%2==0) fscanf(fp,"%d",&L[j*ni2+k]);
	  else fscanf(fp,"%d",&L[(2*j+1)*nind+k]);
	  
	} 
	fscanf(fp,"\n");
      }
    }
    /* printf("Data matrix:\n"); */
    /* showint(L,nind,n2); */
  }
  else{
    if(!(L = (int *) malloc(nind*nloci*sizeof(int)))){ 
      printf("Error in allocating memory in function genotype_s\n\n"); 
      return NULL;
    }
    for(i = 0; i < nind; i++){
      fscanf(fp,"%d",&id[i]); // id nr
      for(j = 0; j < nloci; j++)
	fscanf(fp,"%d",&L[j*nind+i]); 
      fscanf(fp,"\n");
    }
    /* printf("Data matrix:\n"); */
    /* showint(L,nind,nloci); */


  }

  return L;

}
void showint(int *a, int r, int c){
  int i, j;

  for(i = 0; i < r; i++){
    for(j = 0; j < c; j++){

      //cout << fixed << a[j*r+i] << "\t";
      printf(" %d ", a[j*r+i]);
    }
    printf("\n");
  }
}
void show(double *a, int r, int c){
  int i, j;

  for(i = 0; i < r; i++){
    for(j = 0; j < c; j++){

      //cout << fixed << a[j*r+i] << "\t";
      printf("%1.3lf ", a[j*r+i]);
    }
    printf("\n");
  }
}
int allele_count(int *X,int nloci,int nind, int *allcnt,int ploidy,int *Xd){

  int l,i,k,Jl=0,m;
  int *tmp,flag,*ltmp,flag2,m1;
  int ni2=ploidy*nind;
  int nalltot = 0;

  if(!(tmp = (int *) malloc(ploidy*nind*sizeof(int))) ||
     !(ltmp = (int *) malloc(4*sizeof(int)))){ 
    printf("Error in allocating memory in function allele_count\n\n"); 
    return -1;
  }
 
  for(l = 0; l < nloci; l++){ // X is of size nind x (ploidy*nloci)
    ltmp[0] = l*ni2;
    ltmp[1] = (2*l+1)*nind;
    Jl = 1;
    flag2 = 0;
    m1 = 0;
    while(flag2==0){
      if(X[ltmp[0]+m1] != -9 && X[ltmp[0]+m1] != -999){
	tmp[0] = X[ltmp[0]+m1]; // add first allele at loci l
	Xd[ltmp[0]+m1]=0;

	flag2 = 1; // might need to be updated

      }
      else
	Xd[ltmp[0]+m1]=-9;
      
      if(ploidy==2 && X[ltmp[0]+m1]!=X[ltmp[1]+m1] && X[ltmp[1]+m1] != -9 && X[ltmp[1]+m1] != -999){ 
	tmp[1] = X[ltmp[1]+m1];// add second allele at loci l if ind 1 is heterozygote
	Xd[ltmp[1]+m1]=Jl;
	Jl++;
	
      }	 
      else if(ploidy==2 && (X[ltmp[1]+m1] == -9 || X[ltmp[1]+m1] == -999)) Xd[ltmp[1]+m1]=-9;
      m1++;
    }
    /* printf("0.1 add allele %d for ind 1 at loci %d\n",tmp[0],l+1); */
    /* if(Jl==2) printf("0.2 add allele %d for ind 1 at loci %d\n",tmp[1],l+1); */
    for(i = m1; i < nind; i++){
      for(m = 0; m < ploidy; m++){
	flag = 1; // flag for allele l 
	ltmp[m+2] = ltmp[m]+i;
	for(k = 0; k < Jl; k++){ // store allele data for each ind into sparse data matrix Xs
	  /* printf("l = %d, i = %d, k = %d\ttmp[%d] = %d\ttmp[%d] = %d\tX1 = %d\tX2 = %d\n",l,i,k,k,tmp[k],k+nind,tmp[k+nind],X[l*n2+i],X[(2*l+1)*nind+i]); */
	  if(tmp[k] == X[ltmp[m+2]] && (X[ltmp[m+2]] != -9 && X[ltmp[m+2]] != -999)){ 
	    flag = 0;
	    Xd[ltmp[m+2]] = k;
	  }

	  
	  /* printf("Allele %d already existing; allele 1: %d, allele 2: %d\n",k+1,X[l*n2+i],X[(2*l+1)*nind+i]); */

	}
	if((X[ltmp[m+2]] == -9 || X[ltmp[m+2]] == -999)){
	  flag = 0;
	  Xd[ltmp[m+2]] = -9;
	}
	if(flag==1){ // add allele at loci l to temp vector
	  tmp[Jl] = X[ltmp[m+2]];	     
	  /* printf("1. add allele %d for ind %d at loci %d\n",tmp[Jl],i+1,l+1); */
	  Xd[ltmp[m+2]] = Jl;
	  Jl++;
	}
      }
    }
    nalltot += Jl;
    allcnt[l] = Jl;
    // if(Jl==1) ret=0;
    /* printf("number of alleles at loci %d: %d\n",l+1,Jl); */
  }

  free(ltmp);
  free(tmp);
  free(X);

  return nalltot;

}
double *init_alloc(const gsl_rng *gBaseRand, int K,int nind, int nindK, double c, int c0){

  double *Z,*pr;
  int i,k;
  unsigned int *tmp;

  //c = (double)1./K;

  if(!(Z = (double *) malloc(nindK*sizeof(double))) ||
     !(tmp = (unsigned int *) malloc(K*sizeof(int))) ||
     !(pr = (double *) malloc(K*sizeof(double)))){
    printf("Error in allocating memory in function init_alloc\n\n");
    return NULL;
  }

  for(i = 0; i < K; i++) pr[i] = c;
  
  for(i = 0; i < nind; i++){
    /* sample from multinomial distibution to generate allocations for all individuals */
    gsl_ran_multinomial (gBaseRand, K, c0, pr, tmp);
    /* rmultinom(c0, pr, K, tmp); */
    for(k = 0; k < K; k++) Z[k*nind+i] = (double)tmp[k]/c0;
  }

  free(pr);
  free(tmp);

  return Z;

}
int x_ilj(int *Xi, int i, int l, int j, int ploidy, int nind){ /* No of copies of allele j at locus l in ind i */

  int cnt=0,i1;
  
  for(i1 = 0; i1 < ploidy; i1++)
    if(Xi[(ploidy*l+i1)*nind+i]==j) cnt++;
    
  return cnt;

}
int find_max(int *allcnt,int nloci){

  int max=0,i;

  for(i = 0; i < nloci; i++)
    if(allcnt[i]>max)
      max = allcnt[i];
    
  return max;

}
