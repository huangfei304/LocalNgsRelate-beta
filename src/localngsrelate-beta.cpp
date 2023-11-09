#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h> 
#include <libgen.h> 
#include <time.h>
#include <vector>
#include <zlib.h>  
#include <cmath>
#include <cassert>
#include "analysisfunctions.h"
#include "filehandlingfunctions.h"


typedef struct{
  const char *freqfile;
  const char *glafile;
  const char *glbfile;
  const char *outname;
  para p;
}cArg;


void print_info(FILE *fp){
  fprintf(fp, "\n\n");
  fprintf(fp, "You are using LocalNgsRelate version 1.0 with fix a,-k0,-k1,-k2 (build time: %s:%s)\n",__DATE__,__TIME__);
  fprintf(fp, "\nUsage: ./localngsrelate  [options] \n");
  fprintf(fp, "\nRequired options:\n");
  fprintf(fp, "   -f         <filename>   file with frequencies (gz format)\n");
  fprintf(fp, "   -a   <filename>   one file with genotype likelihoods (gz format)\n");
  fprintf(fp,"    -b   <filename>   the other file with genotype likelihoods (gz format)\n");
  fprintf(fp, "   -o         <fileprefix>  Prefix for name of outputfiles\n");
  fprintf(fp, "\n");
  exit(0);
}
double alim[2];
int runoldrelateV;

int main(int argc, char **argv){

  // Start run time clock
  clock_t t=clock();
  time_t t2=time(NULL);
  
  // Set default parameters
  cArg ca;
  ca.freqfile=NULL;
  ca.glafile =NULL;
  ca.glbfile=NULL;
  ca.outname=NULL;
  ca.p.pair[0]=0;
  ca.p.pair[1]=1;
  ca.p.a=-1;
  ca.p.k0=-1;
  ca.p.k1=-1;
  ca.p.k2=-1;

  alim[0]=0.00000001;
  alim[1]=0.10;
  runoldrelateV = 0; // Remove?

  // Specify the expected options for parsing purposes
  static struct option long_options[] = {
      {"f",       required_argument, 0,  'f' },
      {"a",       required_argument, 0, 'a' },
      {"b",       required_argument, 0, 'b' },
      {"o",       required_argument, 0, 'o' },
      {0,        0,                0,  0   }
    };
  
  // Check if any options are specified if not then print calling instruction
  if(argc==1){// if no arguments, print info on program
    print_info(stderr);     
    return 0;
  }
  
  // If there are any options specified then read them
  int opt= 0;
  int long_index = 0;
  while ((opt = getopt_long_only(argc, argv,"", long_options, &long_index )) != -1) {
      switch (opt) {
	// Reading in arguments and setting parameter values accordingly
      case 'f': ca.freqfile = strdup(optarg); break;
      case 'a': ca.glafile = strdup(optarg); break;
      case 'b': ca.glbfile = strdup(optarg); break;
      case 'o': ca.outname = strdup(optarg); break;
      default: return 0;//{fprintf(stderr,"unknown arg:\n");return 0;}
      print_info(stderr);
      }
  }
      
  // Start log file and write command to screen
  if(ca.outname==NULL){
    fprintf(stderr,"\n\t->## Error: the output file is not set\n");
    exit(0);
  }

  std::vector<char *> dumpedFiles;
  fprintf(stderr,"\n\t-> You are using LocalNgsRelate version 1.0 with fix a,-k0,-k1,-k2 (build time: %s:%s)\n",__DATE__,__TIME__);
  fprintf(stderr,"\t-> Command running is:\n\t   ");
  for(int i=0;i<argc;i++){
    fprintf(stderr,"%s ",argv[i]);
  }
  fprintf(stderr,"\n");
  
  // Check if all necessary arguments have been provided
  bool stop = false;
  if(ca.glafile==NULL || ca.glbfile==NULL){
    fprintf(stderr,"\n\t## Error: please supply two genotype likelihood file using the option -a, -b\n");
    stop=true;
  }

  if(ca.freqfile==NULL){
    fprintf(stderr,"\n\t## Error: please supply both a genotype likelihood file and a frequency file\n");
    stop=true;
  }
  if(stop)
    exit(0);

  // Read in freqs
  std::vector<double> freq;
  bgl d = findOverlap(ca.glafile,ca.glbfile, ca.freqfile, freq);
  fprintf(stderr,"\n\t-->It contained frequencies from %d sites for analysis\n",freq.size());

  std::vector<perChr> pd;
  pd = makeDat(d,freq);


  //freq.clear();
  // Print args to a log file
  fprintf(stderr,"\n\t-->mkGenome start\n");
  genome gg = mkGenome(pd,ca.p);

  // If ks are proved by the user and calcA is set: calc a and the corresponding log like
  // a=0.07  k0=0.34 k1=0.35 k2=0.31
  //a=0.10  k0=0.30 k1=0.50 k2=0.20
  double ak012[2][4]={{0.07,0.34,0.35,0.31},{0.10,0.30,0.50,0.20}};
  double *pars;
  double loglike =  calcLoglikeGenome(ak012[0],gg);
  double  curloglike =  calcLoglikeGenome(ak012[1],gg);
  if(loglike>curloglike){
      pars = ak012[0];
  }else{
      pars = ak012[1];
      loglike = curloglike;
  }
  // If one or more parameters need to be estimated this is done

  fprintf(stderr,"\t-> Parameters estimated/set to: a=%.16f k0=%.16f k1=%.16f k2=%.16f loglike=%.16f\n",pars[0],pars[1],pars[2],pars[3],loglike);
  
  // Now do the forward, backward, posterior decoding and viterbi
  if(std::isinf(-loglike)){
    fprintf(stderr,"\t-> Likelihood is 0 with the current parameter values, so no inference was performed.\n");
  }else{
  fprintf(stderr,"\t-> Performing IBD inference along the genome now...\n");
  forward_backward_decode_viterbi(pars,gg);
  
  // Write results to a file
  gzFile add = openFileGz(ca.outname,".IBDmerge.gz","wb",dumpedFiles);
  //gzFile add = openFileGz(ca.outname,".IBDtractinference.gz","wb",dumpedFiles);
  //gzprintf(add,"Chr\tStart\tEnd\tLength\tViterbi\n");
  gzprintf(add,"total\tibd1\tibd1R\tibd2\tibd2R\tkinship\tibd1L\tibd2L\n");
  long int ibdL[5]={0,0,0,0,0};//total,ibd1,ibd2,ibd1L,ibd2L
  for(int i=0;i<pd.size();i++){
    perChr en = pd[i];
    hmmRes to = gg.results[i];
    assert(en.nSites==to.nSites);
    int start= to.pos[0];
    int last = to.pos[0];
    int ibd=to.viterbi[0];
    ibdL[0] += (to.pos[to.nSites-1] - to.pos[0] + 1);
    //fprintf(stderr, "\n\t->chr:%s,start-end:%d-%d,total:%ld\n",en.name,to.pos[0],to.pos[to.nSites-1],ibdL[0]);
    for(int j=1;j<en.nSites;j++){
       if(ibd !=to.viterbi[j]){
         //gzprintf(add,"%s\t%d\t%d\t%d\t%d\n",en.name,start,last,(last-start+1),ibd);
         int len = (last - start + 1 );
         if( ibd > 0 ){
            ibdL[ibd] += len;
			if (len > ibdL[ibd+2])
			  ibdL[ibd+2] = len;
		 }
         start = to.pos[j];
         ibd = to.viterbi[j];
         //gzprintf(add,"%d\t",to.viterbi[j]);
         //gzprintf(add,"%.16e\t%.16e\t%.16e\n",to.post[0][j],to.post[1][j],to.post[2][j]);
       }
       last = to.pos[j];
    }
    //gzprintf(add,"%s\t%d\t%d\t%d\t%d\n",en.name,start,last,(last-start+1),ibd);
	int len = (last - start + 1);
	if(ibd > 0 ){
		ibdL[ibd] +=len;
		if( len > ibdL[ibd+2] )
			ibdL[ibd+2] = len;
	}
  }
  double ibd1R = 1.0*ibdL[1]/ibdL[0]/4.0;
  double ibd2R = 1.0*ibdL[2]/ibdL[0]/2.0;
  double kinship = ibd1R + ibd2R;
  gzprintf(add,"%ld\t%ld\t%f\t%ld\t%f\t%f\t%ld\t%ld\n",ibdL[0],ibdL[1],ibd1R,ibdL[2],ibd2R,kinship,ibdL[3],ibdL[4]);
  gzclose(add);
 }

  fprintf(stderr,"\t-> All analyses are done.\n");
  fprintf(stderr,"\t   [ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr,"\t   [ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  fprintf(stderr,"\n\t   Files created are:\n");
  for(int i=0;1&&i<dumpedFiles.size();i++){
    fprintf(stderr,"\t   %s\n",dumpedFiles[i]);
    free(dumpedFiles[i]);
  }

  
  //delete [] freq ;

  
  return 0;
      
}
