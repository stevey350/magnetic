
#ifndef __MAGNETIC_EDGE_H_
#define __MAGNETIC_EDGE_H_

#include <Eigen/Core>

#include "g2o/core/base_vertex.h"
#include "g2o/core/base_unary_edge.h"

#include "sysdef.h"

namespace Eigen {
typedef Eigen::Matrix<double, 7, 1> Vector7d;
}

// 曲线模型的顶点，模板参数：优化变量维度和数据类型
class VertexParams : public g2o::BaseVertex<7, Eigen::Vector7d >
{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    VertexParams();

    virtual bool read(std::istream& /*is*/);

    virtual bool write(std::ostream& /*os*/) const;

    virtual void setToOriginImpl();                    // 重置

    virtual void oplusImpl(const double* update);      // 更新

};

// 误差模型 模板参数：观测值维度，类型，连接顶点类型
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

