#include <zm/zm_misc.h>

int main(int argc, char *argv[])
{
  double deg, rad;

  if( argc > 1 )
    deg = atof( argv[1] );
  else{
    printf( "angle[deg]> " );
    if( scanf( "%lf", &deg ) == 0 );
  }
  rad = zDeg2Rad( deg );
  printf( "degree %f -> radian %f\n", deg, rad );
  deg = zRad2Deg( rad );
  printf( "radian %f -> degree %f\n", rad, deg );
  return 0;
}
