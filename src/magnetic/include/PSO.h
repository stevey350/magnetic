/*==========================

  Pso.h
  粒子群算法
  应用在初始位置计算

/==========================*/
#ifndef PSO_H
#define PSO_H

#include <stdlib.h>
#include "sysdef.h"

#define object_dim      6
#define max_para_num    7
#define samples         9
#define magnetic_dim    3


//种群个数
const int _pop_num=150;

//迭代次数
const int _max_k=75;

class PSO_Process
{
public:
    PSO_Process();
    ~PSO_Process();

	//初始化       传感器信号      传感器位置    计算初值/结果   传感器个数      目标个数
	void InitData(double (*B)[3], double (*L)[3], double *res, int sensor_num, int object_num)
	{
		_B=B;
		_L=L;
		_res=res;

		_sensor_num=sensor_num;
		_para_num=object_num*object_dim+1;
		_object_num=object_num;
	}

	double Compute();

private:
/*--------------------------

       变量定义

/--------------------------*/
	//传感器信号
	double (*_B)[3];

	//传感器位置
	double (*_L)[3];

	//计算结果
	double *_res;

	//传感器个数
	int _sensor_num;

	//参数个数
	int _para_num;

	//目标个数
	int _object_num;

	//局部最优
	double _best_p[_pop_num][max_para_num+1];

	//全局最优
	double _best_g[max_para_num+1];

	//每一个种群的当前位置
	double _x[_pop_num][max_para_num];

	//每一个种群的当前速度
	double _v[_pop_num][max_para_num];

	//粒子各个参数的取值范围
	double _bound[max_para_num][2];


/*--------------------------

       函数定义

/--------------------------*/

	//产生0~1之间的随机数
	double _rand()
	{
		return (rand()+0.1)/(RAND_MAX+0.1);
	}

	//计算适应度函数
	double _f();

	//
    void SetPara();

};

#endif
