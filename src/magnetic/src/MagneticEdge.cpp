#include <iostream>
#include "MagneticEdge.h"

using namespace std;

VertexParams::VertexParams()
{

}

bool VertexParams::read(std::istream& /*is*/)
{
  cerr << __PRETTY_FUNCTION__ << " not implemented yet" << endl;
  return false;
}

bool VertexParams::write(std::ostream& /*os*/) const
{
  cerr << __PRETTY_FUNCTION__ << " not implemented yet" << endl;
  return false;
}

void VertexParams::setToOriginImpl()                    // 重置
{
  _estimate << 0, 0, 0, 0, 0, 0, 0;
}

void VertexParams::oplusImpl(const double* update)      // 更新
{
  Eigen::Vector7d::ConstMapType v(update);
  _estimate += v;

  // constraints for m^2+n^2+p^2 = 1;
//      double one = sqrt(_estimate(3)*_estimate(3) +
//                        _estimate(4)*_estimate(4) +
//                        _estimate(5)*_estimate(5));
//      _estimate(3) = _estimate(3)/one;
//      _estimate(4) = _estimate(4)/one;
//      _estimate(5) = _estimate(5)/one;
}


//-----------------------------------------------------------------------------------------------------
EdgePointOnCurve::EdgePointOnCurve( Eigen::Vector3d position): BaseUnaryEdge(), _sensor_position(position)
{

}

bool EdgePointOnCurve::read(std::istream& /*is*/)
{
  cerr << __PRETTY_FUNCTION__ << " not implemented yet" << endl;
  return false;
}

bool EdgePointOnCurve::write(std::ostream& /*os*/) const
{
  cerr << __PRETTY_FUNCTION__ << " not implemented yet" << endl;
  return false;
}

void EdgePointOnCurve::computeError()
{
  // estimated value (a, b, c, m, n, p, Bt)
  const VertexParams* params = static_cast<const VertexParams*>(vertex(0));
  const double& a  = params->estimate()(0);
  const double& b  = params->estimate()(1);
  const double& c  = params->estimate()(2);
  const double& m  = params->estimate()(3);
  const double& n  = params->estimate()(4);
  const double& p  = params->estimate()(5);
  const double& Bt = params->estimate()(6);

  // position of sensor
  const double &x = _sensor_position(0);
  const double &y = _sensor_position(1);
  const double &z = _sensor_position(2);

  // function value based on dipole model
  double R = sqrt( (x-a)*(x-a) + (y-b)*(y-b) + (z-c)*(z-c) );
  double fval0 = Bt*((3*(m*(x-a)+n*(y-b)+p*(z-c))*(x-a))/(R*R*R*R*R) - m/(R*R*R));
  double fval1 = Bt*((3*(m*(x-a)+n*(y-b)+p*(z-c))*(y-b))/(R*R*R*R*R) - n/(R*R*R));
  double fval2 = Bt*((3*(m*(x-a)+n*(y-b)+p*(z-c))*(z-c))/(R*R*R*R*R) - p/(R*R*R));

  _error(0) = 0.5*(fval0 - measurement()(0))*(fval0 - measurement()(0));
  _error(1) = 0.5*(fval1 - measurement()(1))*(fval1 - measurement()(1));
  _error(2) = 0.5*(fval2 - measurement()(2))*(fval2 - measurement()(2));
}
/*
void EdgePointOnCurve::linearizeOplus()
{
    // estimated value (a, b, c, m, n, p)
    const VertexParams* params = static_cast<const VertexParams*>(vertex(0));
    const double& a = params->estimate()(0);
    const double& b = params->estimate()(1);
    const double& c = params->estimate()(2);
    const double& m = params->estimate()(3);
    const double& n = params->estimate()(4);
    const double& p = params->estimate()(5);
    const double& Bt = params->estimate()(6);

    // position of sensor
    const double &x = _sensor_position(0);
    const double &y = _sensor_position(1);
    const double &z = _sensor_position(2);

    double R = sqrt((x-a)*(x-a)+(y-b)*(y-b)+(z-c)*(z-c));
    double xa, yb, zc;
    xa = x-a; yb = y-b; zc = z-c;
    double R3, R5, R7;
    R3=R*R*R; R5=R3*R*R; R7=R5*R*R;

    double tmp1=3*m*(x-a)+3*n*(y-b)+3*p*(z-c);

    // Bx对a, b, c, m, n, p求偏导
    _jacobianOplusXi(0,0) = Bt*((-3*m*xa-tmp1-3*m*xa)/R5+5*tmp1*xa*xa/R7);
    _jacobianOplusXi(0,1) = Bt*((-3*n*xa     -3*m*yb)/R5+5*tmp1*xa*yb/R7);
    _jacobianOplusXi(0,2) = Bt*((-3*p*xa     -3*m*zc)/R5+5*tmp1*xa*zc/R7);
    _jacobianOplusXi(0,3) = Bt*(3*xa*xa/R5-1/R3);
    _jacobianOplusXi(0,4) = Bt* 3*yb*xa/R5;
    _jacobianOplusXi(0,5) = Bt* 3*zc*xa/R5;
    _jacobianOplusXi(0,6) = tmp1*xa/R5-m/R3;

    // By对a, b, c, m, n, p求偏导
    _jacobianOplusXi(1,0) = Bt*((-3*m*yb     -3*n*xa)/R5+5*tmp1*xa*yb/R7);
    _jacobianOplusXi(1,1) = Bt*((-3*n*yb-tmp1-3*n*yb)/R5+5*tmp1*yb*yb/R7);
    _jacobianOplusXi(1,2) = Bt*((-3*p*yb     -3*n*zc)/R5+5*tmp1*zc*yb/R7);
    _jacobianOplusXi(1,3) = Bt* 3*yb*xa/R5;
    _jacobianOplusXi(1,4) = Bt*(3*yb*yb/R5-1/R3);
    _jacobianOplusXi(1,5) = Bt* 3*yb*zc/R5;
    _jacobianOplusXi(1,6) = tmp1*yb/R5-n/R3;

    // Bz对a, b, c, m, n, p求偏导
    _jacobianOplusXi(2,0) = Bt*((-3*m*zc     -3*p*xa)/R5+5*tmp1*zc*xa/R7);
    _jacobianOplusXi(2,1) = Bt*((-3*n*zc     -3*p*yb)/R5+5*tmp1*zc*yb/R7);
    _jacobianOplusXi(2,2) = Bt*((-3*p*zc-tmp1-3*p*zc)/R5+5*tmp1*zc*zc/R7);
    _jacobianOplusXi(2,3) = Bt* 3*zc*xa/R5;
    _jacobianOplusXi(2,4) = Bt* 3*zc*yb/R5;
    _jacobianOplusXi(2,5) = Bt*(3*zc*zc/R5-1/R3);
    _jacobianOplusXi(2,6) = tmp1*zc/R5-p/R3;
}
*/


EdgeOrientation::EdgeOrientation(): BaseUnaryEdge()
{

}

bool EdgeOrientation::read(std::istream& /*is*/)
{
  cerr << __PRETTY_FUNCTION__ << " not implemented yet" << endl;
  return false;
}

bool EdgeOrientation::write(std::ostream& /*os*/) const
{
  cerr << __PRETTY_FUNCTION__ << " not implemented yet" << endl;
  return false;
}

#define D            0.076
#define r_c          0.005           // 5 mm
#define r_p          0.1
void EdgeOrientation::computeError()
{
  // estimated value (a, b, c, m, n, p, Bt)
  const VertexParams* params = static_cast<const VertexParams*>(vertex(0));
  const double& c  = params->estimate()(2);
  const double& m  = params->estimate()(3);
  const double& n  = params->estimate()(4);
  const double& p  = params->estimate()(5);
  const double& Bt = params->estimate()(6);

//  double max_value = (Bt>1e-9 ? Bt : 1e-9);
  double min_value = Bt*1e8>0 ? 0 : Bt*1e8;
  double max_c = abs(c-D) > r_c ? abs(c-D) : r_c;
  double max_p = abs(p) > r_p ? abs(p) : r_p;
  float sigma1, sigma2;

  sigma1 = 3e-3;                            // Before improvement -> sigma1 = 3e-6;
  sigma2 = 2e-3;                            // Before improvement -> sigma2 = 2e-9


  _error(0) = 0.5*(m*m + n*n + p*p - 1)*(m*m + n*n + p*p - 1);
//  _error(1) = Bt-7.918e-8;            // BT is regarded as a constant
//  _error(1) = 1e-17/(max_value*fabs(Bt)*1e8);

//  _error(1) = 0.5*min_value*min_value;
//  _error(2) = Bt<1e-8 ? 1e-13/Bt: 0;

//  _error(3) = sigma1*0.5*(c-0.076)*(c-0.076);
//  _error(4) = sigma2*0.5*((p-0)*(p-0));

//  _error(3) = sigma1*0.5*(max_c-r_c)*(max_c-r_c);           // (c-0.076)*(c-0.076);
//  _error(4) = sigma2*0.5*(max_p-r_p)*(max_p-r_p);           //   sigma2*0.5*(p-0)*(p-0);

}

