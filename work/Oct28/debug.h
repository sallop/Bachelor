#ifndef __INCLUDE_DEBUG_H__
#define __INCLUDE_DEBUG_H__

#include <cstdarg>
#include <cstring>

extern int sat[1+6], at[1+6];
extern int nnn, j, nt;
extern double t;

FILE *fp_debug;

#define dbg(fmt, ...) debug_print(__FILE__, __FUNCTION__, __LINE__, (fmt), __VA_ARGS__)

void init_debug(const char *ofname)
{
  fp_debug = fopen(ofname, "w");
}

void finalize_debug()
{
  fclose(fp_debug);
}

void debug_print(const char *file,
		 const char *func,
		 int line,
		 const char *fmt,
		 ...)
{
  char buf1[32], buf2[32];
  va_list args;
  va_start(args, fmt);
  vsprintf(buf1, fmt, args);
  va_end(args);

  sprintf(buf2, "%s:%s:%d:%s\n", file, func, line, buf1);

  fprintf(fp_debug, buf2);
}

FILE *fp_assert_tansaku;
FILE *fp_assert_time;
FILE *fp_assert_repeat;

void init_assertion()
{
  fp_assert_tansaku = fopen("DEBUG/tansaku.dat","w");
  fp_assert_time = fopen("DEBUG/time.dat","w");
  fp_assert_repeat = fopen("DEBUG/repeat.dat","w");
}

//void assert_tansaku()
//{
//  fprintf(fp_assert_tansaku, "nnn = %d, j = %d\n", nnn, j);
//}

//void assert_time(int flag)
//{
//  if (flag)
//    fprintf(fp_assert_time,"nt = %d, t = %lf\n", nt, t);
//}

void assert_repeat()
{
  fprintf(fp_assert_repeat, "nnn=%d,", nnn);
  fprintf(fp_assert_repeat, "sat[1]=%d,", sat[1]);
  fprintf(fp_assert_repeat, "sat[2]=%d,", sat[2]);
  fprintf(fp_assert_repeat, "sat[3]=%d,", sat[3]);
  fprintf(fp_assert_repeat, "sat[4]=%d\n", sat[4]);
}

void finalize_assertion()
{
  fclose(fp_assert_tansaku);
  fclose(fp_assert_time);
  fclose(fp_assert_repeat);
}

#endif // __INCLUDE_DEBUG_H__
