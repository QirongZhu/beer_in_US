namespace Kepler{

  const double UNIVERSAL_PI = atan(1.0) * 4;
  const int SUCCESS = 1;
  const int FAILURE = 0;
  
  struct State
  {
    double x, y, z, xd, yd, zd;

    State(): x(0.0), y(0.0), z(0.0), xd(0.0), yd(0.0), zd(0.0){}

    State(double x_, double y_, double z_,
	  double vx_, double vy_, double vz_):
      x(x_), y(y_), z(z_), xd(vx_), yd(vy_), zd(vz_)  {}

    friend void kepler_step(State &s, double dt);
    friend void tidal_kick(State &s, double dt);
  };  
  
}
