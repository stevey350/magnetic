
#ifndef __MAGNETIC_EDGE_H_
#define __MAGNETIC_EDGE_H_

#include <Eigen/Core>

#include "g2o/core/base_vertex.h"
#include "g2o/core/base_unary_edge.h"

#include "sysdef.h"

namespace Eigen {
typedef Eigen::Matrix<double, 7, 1> Vector7d;
}

// ����ģ�͵Ķ��㣬ģ��������Ż�����ά�Ⱥ���������
class VertexParams : public g2o::BaseVertex<7, Eigen::Vector7d >
{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    VertexParams();

    virtual bool read(std::istream& /*is*/);

    virtual bool write(std::ostream& /*os*/) const;

    virtual void setToOriginImpl();                    // ����

    virtual void oplusImpl(const double* update);      // ����

};

// ���ģ�� ģ��������۲�ֵά�ȣ����ͣ����Ӷ�������
class EdgePointOnCurve : public g2o::BaseUnaryEdge<SENSOR_DIM, Eigen::Vector3d, VertexParams>
{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    EdgePointOnCurve( Eigen::Vector3d position);    // : BaseUnaryEdge(), _sensor_position(position)
    virtual bool read(std::istream& /*is*/);
    virtual bool write(std::ostream& /*os*/) const;

    void computeError();

//    void linearizeOplus();

private:
    Eigen::Vector3d _sensor_position;
};

class EdgeOrientation : public g2o::BaseUnaryEdge<5, Eigen::Vector3d, VertexParams>
{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    EdgeOrientation();  // : BaseUnaryEdge()
    virtual bool read(std::istream& /*is*/);
    virtual bool write(std::ostream& /*os*/) const;

    void computeError();
};

#endif

