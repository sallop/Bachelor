#include <cstdio>
#include <unistd.h>
#include <stdlib.h>
#include <iostream>

enum SeeSide{
  normal, xside, zside,
};

void
Animation(SeeSide see, int start, int end)
{
  FILE *gp;
  int i, k;
  gp = popen("gnuplot","w");

  fprintf(gp,"set xrange[-0.15:0.15]\n");
  fprintf(gp,"set yrange[-0.07:0.20]\n");
  fprintf(gp,"set zrange[ 0.00:0.30]\n");  
  fprintf(gp,"set terminal x11\n");

  switch(see){
  case normal:
    fprintf(gp, "set view 45, 300\n");
    break;
  case xside:
    fprintf(gp, "set view 90, 0\n");
    break;
  case zside:
    fprintf(gp, "set view 0, 0\n");
    break;
  default:
    exit(EXIT_FAILURE);
  }
  fprintf(gp,"set size 0.6, 0.6\n");
  fprintf(gp," set grid \n");

  for(k = 1;k <= 2; k++){
    for(i = start; i <= end; i += 10){    
      fprintf(gp,
	      "splot "
	      "'dir-dat/s1-%05d.dat' w l 1,"
	      "'dir-dat/s2-%05d.dat' w l 3,"
	      "'dir-dat/s3-%05d.dat' w l 4,"
	      "0.10 5\n"
	      , i, i, i);
	      //"'dir-dat/e1-%05d.dat' w l 3\n", i, i);
      //"\n",i,i);
      fprintf(gp,"pause 0.05 \n");
    }
  }
  //  for(i = 10; i <= 640; i += 10){
  for(i = start; i <= end; i += 10){
    fprintf(gp,
	    "splot "
	    "'dir-dat/s1-%05d.dat' w l 1,"
	    "'dir-dat/s2-%05d.dat' w l 3,"
	    "'dir-dat/s3-%05d.dat' w l 4,"
	    "0.10 5\n"
	    ,i,i,i);
	    //"'dir-dat/e1-%05d.dat' w l 3\n",i,i);
    //x+0.05 4,-x+0.05 3\n",i,i);
    fprintf(gp,"set terminal gif\n");
    fprintf(gp,"set grid\n");
    fprintf(gp,"set output 'dir-gif/work-2_%05d.gif'\n",i);
    fprintf(gp,"rep\n");
  }
  
  for(k = 1;k <= 2; k++){
    for(i = start; i <= end; i += 10){      
      fprintf(gp,
	      "splot "
	      "'dir-dat/s1-%05d.dat' w l 1,"
	      "'dir-dat/s2-%05d.dat' w l 3,"
	      "'dir-dat/s3-%05d.dat' w l 4,"
	      "0.10 5\n"
	      ,i,i,i);	      
      //      "'dir-dat/e2-%05d.dat' w l 3\n"
      //	      , i, i);
      //"\n",i,i);
      fprintf(gp,"pause 0.05 \n");
    }
  }
  for(i = start; i <= end; i += 10){
    fprintf(gp,
	    "splot "
	    "'dir-dat/s1-%05d.dat' w l 1,"
	    "'dir-dat/s2-%05d.dat' w l 3,"
	    "'dir-dat/s3-%05d.dat' w l 4,"
	    "0.10 5\n"
	    ,i,i,i);
    //"'dir-dat/e2-%05d.dat' w l 3\n",i,i);
    //x+0.05 4,-x+0.05 3\n",i,i);
    fprintf(gp,"set terminal gif\n");
    fprintf(gp,"set grid\n");
    fprintf(gp,"set output 'dir-gif/work-2.gif'\n",i);
    fprintf(gp,"rep\n");
  }

  pclose(gp);  
}

int
main(int argc, char **argv)
{
  int opt;
  int start, end, time;
  enum SeeSide flag;

  start = 0; end = 640;
  while((opt = getopt(argc, argv, "s:e:nxz")) != -1){
    switch(opt){
    case 's':
      start = atoi(optarg);
      break;
    case 'e':
      end = atoi(optarg);
      break;      
    case 'n':
      flag = normal;
      break;
    case 'x':
      flag = xside;
      break;
    case 'z':
      flag = zside;
      break;
    default:			// '?'
      fprintf(stderr, "Usage: %s\n", argv[0]);
      exit(EXIT_FAILURE);
    }
  }

  printf("start=%d,end=%d\n", start, end);
  Animation(flag, start, end);

  return 0;
}
