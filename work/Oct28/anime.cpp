#include <cstdio>
#include <cstring>
#include <cmath>

int main(int argc, char *argv[])
{
  char fname[32] = "";
  
  if(!strcmp(argv[1],"dp")){
    sprintf(fname,"dir-dat/dp");
  }

  if(!strcmp(argv[1],"e1")){
    sprintf(fname,"dir-dat/e1");
  }  

  FILE *gp = popen("gnuplot","w");
  fprintf(gp, "set xrange[-0.2:0.2]\n");
  fprintf(gp, "set yrange[-0.2:0.2]\n");
  fprintf(gp, "set zrange[ 0.0:0.4]\n");              
  for(int i = 0; i < 630; i+=10){
    char f[32];
    sprintf(f, fname, i);
    fprintf(gp, "splot '%s-%05d.dat' w l\n", f, i);
  }
  fclose(gp);

  
  return 0;
}
