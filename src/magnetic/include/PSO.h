/*==========================

  Pso.h
  ����Ⱥ�㷨
  Ӧ���ڳ�ʼλ�ü���

/==========================*/
#ifndef PSO_H
#define PSO_H

#include <stdlib.h>
#include "sysdef.h"

#define object_dim      6
#define max_para_num    7
#define samples         9
#define magnetic_dim    3


//��Ⱥ����
const int _pop_num=150;

//��������
const int _max_k=75;

class PSO_Process
{
public:
    PSO_Process();
    ~PSO_Process();

	//��ʼ��       �������ź�      ������λ��    �����ֵ/���   ����������      Ŀ�����
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

       ��������

/--------------------------*/
	//�������ź�
	double (*_B)[3];

	//������λ��
	double (*_L)[3];

	//������
	double *_res;

	//����������
	int _sensor_num;

	//��������
	int _para_num;

	//Ŀ�����
	int _object_num;

	//�ֲ�����
	double _best_p[_pop_num][max_para_num+1];

	//ȫ������
	double _best_g[max_para_num+1];

	//ÿһ����Ⱥ�ĵ�ǰλ��
	double _x[_pop_num][max_para_num];

	//ÿһ����Ⱥ�ĵ�ǰ�ٶ�
	double _v[_pop_num][max_para_num];

	//���Ӹ���������ȡֵ��Χ
	double _bound[max_para_num][2];


/*--------------------------

       ��������

/--------------------------*/

	//����0~1֮��������
	double _rand()
	{
		return (rand()+0.1)/(RAND_MAX+0.1);
	}

	//������Ӧ�Ⱥ���
	double _f();

	//
    void SetPara();

};

#endif
