// g2o - General Graph Optimization
// Copyright (C) 2012 R. Kümmerle
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "ros/ros.h"
#include "std_msgs/String.h"

#include <Eigen/Core>
#include <math.h>
#include <iostream>

#include "SensorsData.h"
#include "AsyncSerial.h"
#include "sysdef.h"
#include "MagneticEdge.h"
#include "key_board.h"

#include "g2o/stuff/sampler.h"
#include "g2o/stuff/command_args.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/solver.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/core/base_vertex.h"
#include "g2o/core/base_unary_edge.h"
#include "g2o/solvers/dense/linear_solver_dense.h"

#include "LM_Process.h"

using namespace std;

int main(int argc, char** argv)
{
    // 1-ros initial
    ros::init(argc, argv, "magnetic_node");
    ros::NodeHandle nh;

    ros::Publisher magnetic_pub = nh.advertise<std_msgs::String>("magnetic_pub", 10);
    ros::Rate loop_rate(20);

    // 2-serial initial for reading sensors data
    CallbackAsyncSerial serial;

    try
    {
        serial.open("/dev/ttyACM0",115200);
    }
    catch(boost::system::system_error& e)
    {
        //Errors during open
        cout<<"Error: "<<e.what()<< " Please check the connection of STM32 module." << endl;
        return 1;
    }

    serial.setCallback(bind(readCallback, _1, _2));//&

    // 3-key board initial
    tty_set();
    int save_flag = 0;


    bool verbose = false;
    double w_sigma = 0.02;                 // 噪声Sigma值
    int maxIterations = 300;

    Eigen::Matrix<double, SENSOR_NUM, SENSOR_DIM> Bm;
    //0.000000, 0.000000, 0.106000, 0.000000, 0.000000, 1.000000.
    Bm <<  -0.00002973, 0.00002972, 0.00001833,
            -0.00000000, 0.00005621, 0.00004664,
            0.00002971, 0.00002957, 0.00001818,
            -0.00005620, -0.00000000, 0.00004661,
            -0.00000000, -0.00000000, 0.00013805,
            0.00005608, -0.00000000, 0.00004619,
            -0.00002975, -0.00002972, 0.00001834,
            -0.00000000, -0.00005622, 0.00004668,
            0.00002973, -0.00002957, 0.00001819;

    // some handy typedefs
    typedef g2o::BlockSolver< g2o::BlockSolverTraits<Eigen::Dynamic, Eigen::Dynamic> >  MyBlockSolver;
    typedef g2o::LinearSolverDense<MyBlockSolver::PoseMatrixType> MyLinearSolver;

    auto linearSolver = g2o::make_unique<MyLinearSolver>();
    auto blockSolver = g2o::make_unique<MyBlockSolver>(std::move(linearSolver));

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(std::move(blockSolver));

    // setup the solver
    g2o::SparseOptimizer optimizer;

    optimizer.setAlgorithm(solver);

    // 1. add the parameter vertex
    Eigen::Vector7d init_value;
    init_value << -0.03, -0.03, 0.076, 0, 0, 1, Bt_guess;
    VertexParams* params = new VertexParams();
    params->setId(0);
    params->setEstimate(init_value); // some initial value for the params
    optimizer.addVertex(params);

    // 2. add the points we measured to be on the curve
    EdgePointOnCurve* edge_point[SENSOR_NUM];
    for (int i = 0; i < SENSOR_NUM; ++i)
    {
        // sensor position
        Eigen::Vector3d position;
        position << sensor_position[i][0], sensor_position[i][1], sensor_position[i][2];

        // sensor measured values
        Eigen::Vector3d Bl;
        Bl << Bm(i, 0), Bm(i, 1), Bm(i, 2);
     //   Bl << sensors_data[i][0] * 1e-6, sensors_data[i][1] * 1e-6, sensors_data[i][2] * 1e-6;

        edge_point[i] = new EdgePointOnCurve(position);
        edge_point[i]->setInformation(Eigen::Matrix3d::Identity() );     //xxx
        edge_point[i]->setVertex(0, params);
        edge_point[i]->setMeasurement(Bl);
        optimizer.addEdge(edge_point[i]);
    }

    EdgeOrientation *edge = new EdgeOrientation();
    edge->setInformation(Eigen::Matrix<double, 5, 5>::Identity() );
    edge->setVertex(0, params);
    optimizer.addEdge(edge);

    float a, b, c, m, n, p, bt;

    // LM
//    LM_Process obj_res;
//    double object_para[7] = {0.003, 0.003, 0.076, 0, 0.1, 0.9, Bt_guess};

    while(ros::ok())
    {
        if(0 == read_frame_ok)
        {
            usleep(1000);
            continue;
        }
        read_frame_ok = 0;

        if(0 == geomagnet_removed_ok)
            continue;

        double Bst = 0;
        for (int i=0; i<SENSOR_NUM; i++)
        {
            Bst += fabs(valid_sensors_data[i][0]) + fabs(valid_sensors_data[i][1]) + fabs(valid_sensors_data[i][2]);
        }
//        cout << "Bst = " << Bst << endl;
        if(Bst <= 670e-6)
            continue;

        /*
        obj_res.InitData(valid_sensors_data, sensor_position, object_para, 9, 1);
        obj_res.Compute();

        a = object_para[0];
        b = object_para[1];
        c = object_para[2];
        m = object_para[3];
        n = object_para[4];
        p = object_para[5];
        bt = object_para[6];
        cout << "abcmnp=" << a << " " << b << " " << c << " " << m << " " << n << " " << p << " " << bt << endl;
        */


        // set measurement values
        for (int i=0; i<SENSOR_NUM; i++)
        {
            Eigen::Vector3d Bl;
         //   Bl << Bm(i, 0), Bm(i, 1), Bm(i, 2);
            Bl << valid_sensors_data[i][0], valid_sensors_data[i][1], valid_sensors_data[i][2];

            edge_point[i]->setMeasurement(Bl);
        }

        // 执行优化
//        cout<<"start optimization"<<endl;
        chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
        optimizer.initializeOptimization();
        optimizer.setVerbose(verbose);
        optimizer.optimize(maxIterations);
        chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
        chrono::duration<double> time_used = chrono::duration_cast<chrono::duration<double>>( t2-t1 );
//        cout<<"solve time cost = "<<time_used.count()<<" seconds. "<<endl;

        if (verbose)
            cout << endl;

        // print out the result
//        cout << "Iterative least squares solution" << endl;
        a = params->estimate()(0);
        b = params->estimate()(1);
        c = params->estimate()(2);
        m = params->estimate()(3);
        n = params->estimate()(4);
        p = params->estimate()(5);
        bt = params->estimate()(6);
        cout << "abcmnp=" << a << " " << b << " " << c << " " << m << " " << n << " " << p << " " << bt << endl;
        cout << "theta = " << atan2(n, m)*180/M_PI << endl;
//        cout << "abcmnp = " << params->estimate()(0) << " " << params->estimate()(1) << " " << params->estimate()(2) << " " <<
//             params->estimate()(3) << " " << params->estimate()(4) << " " << params->estimate()(5) << " " << params->estimate()(6) << endl;

        // set initial values
        Eigen::Vector7d next_value;
        next_value << a, b, c, m, n, p, bt;
        params->setEstimate(next_value);

        // scan key board
        if( kbhit() )
        {
            const int key = getchar();
            if(key == 's')
            {
                save_flag = 10;
                cout << "start save data." << endl;
            }

            if(key == 'q')
                break;
        }

        // save data
        if(save_flag > 0)
        {
            fstream outfile;
            outfile.open("data.txt", ios::binary | ios::app | ios::in | ios::out);
            outfile.precision(5);

            outfile << a << " " << b << " " << c << " " << m << " " << n << " " << p << " " << bt << endl;

            outfile.close();
            save_flag--;

            if(0 == save_flag)
                cout << "save OK!" << endl;
        }

        //
        std_msgs::String msg;

        std::stringstream ss;
        ss << "value=" << next_value;

        msg.data = ss.str();

        magnetic_pub.publish(msg);

        ros::spinOnce();

        loop_rate.sleep();
    }

    tty_reset();

    return 0;
}
