
/*********************

LM_Process.h

�������޸�Ϊ2ά����

************************/
#ifndef LM_PROCESS_H
#define LM_PROCESS_H

#include "sysdef.h"

#define object_dim      6
#define max_para_num    7
#define samples         9
#define magnetic_dim    3


class LM_Process
{
public:

    LM_Process();
    ~LM_Process();

	void InitData(double (*B)[3], double (*L)[3], double* res, int sensor_num, int object_num)
	{
		_B=B;
		_L=L;
		_res=res;

		_sensor_num=sensor_num;
		_object_num=object_num;
		_para_num=object_num*object_dim+1;
	};
	double Compute();

private:

	//===============================================
	//��������   ���б������� _ ��ͷ
    //===============================================

	//�������ź�
	double (*_B)[3];
	
	//����������
	double (*_L)[3];

	//�����ʼֵ�Լ����
	double* _res;

	//   ����������   ��������   Ŀ�����
	int _sensor_num, _para_num, _object_num;

	//���������дų�����ֵ B_guess=f(res)
	//double _B_guess[samples * magnetic_dim];

	//��������ֵ res_new=res+Ep
	double _res_new[max_para_num];

	//�ų�ģ�ͺ��� f ���ſɱȾ���
	double _J[samples * magnetic_dim][max_para_num]; 
	
	//A=JT*J
	double _A[max_para_num][max_para_num];

	//A��copy
	double _A2[max_para_num][max_para_num];

	//Ep=B-B_guess  
	double _Ep[samples * magnetic_dim];
	double _Ep2[samples * magnetic_dim];

	//Sp=G/(A+uI) 
	double _Sp[max_para_num];

	//G=JT*Ep
	double _G[max_para_num];

	//���ƶ�����Ϣ
	double _u[max_para_num];


	//==============================================
	//��������
	//==============================================

	void set_Ep(double *Ep, double *res);
	void set_J();

	//calculate A and G
	void set_A();

	//calculate Sp
	bool solve();

	// ||vlaue||
	double model(double *value, int len); 

	// u�ĵ���ֵ
	double adjust();


};

#endif
