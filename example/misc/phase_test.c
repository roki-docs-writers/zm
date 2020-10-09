#include <zm/zm_misc.h>

#define DT 0.001

int main(void)
{
  double p1, p2;
  register int i;

  for( i=-5000; i<=5000; i++ ){
    p1 = 2*zPI*DT*i;
    p2 = zPhaseNormalize( p1 );
    printf( "%f %f\n", p1, p2 );
  }
  return 0;
}
