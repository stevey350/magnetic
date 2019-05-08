
/*=========================

  PSO_Process.cpp

/=========================*/
#include<math.h>
#include "PSO.h"


PSO_Process::PSO_Process()
{

}

PSO_Process::~PSO_Process()
{

}

/*==================================

  成员函数SetPara
  设置bound值
  设置x，v的初始值

/==================================*/
void PSO_Process::SetPara()
{
	int i,j,k;

	//设置bound
	//for (i=0;i<_para_num-1;i++)
	for (i=0;i<_object_num;i++)
	{
		/*j=i%6;
		
		_bound[i][0]=j>2 ? -1 : -0.4;
		_bound[i][1]=j>2 ?  1 :  0.4;*/

		j=i*object_dim;

		_bound[0+j][0]=-0.3; _bound[0+j][1]= 0.3;
		_bound[1+j][0]=-0.3; _bound[1+j][1]= 0.3;
		_bound[2+j][0]= 0;   _bound[2+j][1]= 0.3;
		_bound[3+j][0]=-1  ; _bound[3+j][1]= 1  ;
		_bound[4+j][0]=-1  ; _bound[4+j][1]= 1  ;
		_bound[5+j][0]=-1  ; _bound[5+j][1]= 1  ;
	}
    _bound[_para_num-1][0]=Bt_guess-0.3e-7; // 原来是bt_guess-0.1
    _bound[_para_num-1][1]=Bt_guess+0.3e-7;	     // 原来是bt_guess+0.1

	_best_g[_para_num]=3e+15;

	//设置x，v
	for (i=0;i<_pop_num;i++)
	{
		for (j=0;j<_para_num;j++)
		{
			_v[i][j]=_rand()*(_bound[j][1]-_bound[j][0])/2;
			_x[i][j]=_rand()*(_bound[j][1]-_bound[j][0])+_bound[j][0];

			_best_p[i][j]=_x[i][j];
			_res[j]=_x[i][j];
		}
		_best_p[i][_para_num]=_f();
		
		if(_best_g[_para_num]>_best_p[i][_para_num])
		{
			_best_g[_para_num]=_best_p[i][_para_num];
			k=i;
		}
	}

	for (i=0;i<_para_num;i++)
		_best_g[i]=_best_p[k][i];

}

/*==================================

  成员函数f
  返回适应度

/==================================*/
double PSO_Process::_f()
{
	double res=0;
	double a,b,c,m,n,p,bt,x,y,z,R;
	int i,j,k;
	double bx,by,bz;

	bt=_res[ object_dim * _object_num ];

	for (i=0;i<_sensor_num;i++)
	{
		x=_L[i][0];
		y=_L[i][1];
		z=_L[i][2];

		bx=_B[i][0];
		by=_B[i][1];
		bz=_B[i][2];

		for (j=0;j<_object_num;j++)
		{
        
    		k=object_dim*j;

    		a=_res[0+k];b=_res[1+k];c=_res[2+k];
        	m=_res[3+k];n=_res[4+k];p=_res[5+k];

    		R=sqrt((x-a)*(x-a)+(y-b)*(y-b)+(z-c)*(z-c));

    		bx-=bt*((3*(m*(x-a)+n*(y-b)+p*(z-c))*(x-a))/(R*R*R*R*R)-m/(R*R*R));
    		by-=bt*((3*(m*(x-a)+n*(y-b)+p*(z-c))*(y-b))/(R*R*R*R*R)-n/(R*R*R));
    		bz-=bt*((3*(m*(x-a)+n*(y-b)+p*(z-c))*(z-c))/(R*R*R*R*R)-p/(R*R*R));

		}

		res+=bx*bx+by*by+bz*bz;
	}

	return sqrt(res);
}

/*==================================

  成员函数Compute
  迭代计算目标参数 

/==================================*/
double PSO_Process::Compute()
{
	int i,j,k;
	int count=0;
	double w=2;
	double c1=2,c2=2;
	bool stop=false;
	bool change_best=false;
	double tmp;
    int outflag;

	SetPara();

	while(!stop)
	{
		w*=0.8;
		for(i=0;i<_pop_num;i++)
		{
			for (j=0;j<_para_num;j++)
			{				
				_v[i][j]=w*_v[i][j]+c1*_rand()*(_best_p[i][j]-_x[i][j])+c2*_rand()*(_best_g[j]-_x[i][j]);
				//_v[i][j]=_v[i][j]>_bound[i][1] ? _bound[i][1] : _v[i][j];
				//_v[i][j]=_v[i][j]<_bound[i][0] ? _bound[i][0] : _v[i][j];

				_x[i][j]=_x[i][j]+_v[i][j];
				outflag= (_x[i][j]<_bound[j][0] || _x[i][j]>_bound[j][1]);
				_x[i][j]=_x[i][j]-outflag*_v[i][j];
				
			/*	if (_x[i][j]>_bound[j][1])
					_x[i][j]=_bound[j][1]-0.1;

                if (_x[i][j]<_bound[j][0])
					_x[i][j]=_bound[j][0]+0.1; */

			}//for (j=0;j<_para_num;j++)

			for (j=0;j<_para_num;j++)
				_res[j]=_x[i][j];

			for (j=0;j<_object_num;j++)
			{
				double ones=sqrt(_res[3+j*6]*_res[3+j*6]+_res[4+j*6]*_res[4+j*6]+_res[5+j*6]*_res[5+j*6]);
				_res[3+j*6]/=ones;
				_res[4+j*6]/=ones;
				_res[5+j*6]/=ones;
			}

			tmp=_f();
			if (tmp<_best_p[i][_para_num])
			{
				_best_p[i][_para_num]=tmp;
				for (j=0;j<_para_num;j++)
					_best_p[i][j]=_x[i][j];
			}

			if(_best_p[i][_para_num]<_best_g[_para_num])
			{
				k=i;
				_best_g[_para_num]=_best_p[i][_para_num];
				change_best=true;
			}
		}//for(i=0;i<_pop_num;i++)

		if(change_best)
		{
			for (i=0;i<_para_num;i++)
				_best_g[i]=_best_p[k][i];
			change_best=false;
		}

		//stop=(fabs(_best_g[_para_num])<1000 || ++count>_max_k) ? true : false ;
		stop=(fabs(_best_g[_para_num])<1.0e-5 || ++count>_max_k) ? true : false ;

	}//while(!stop)

	for(i=0;i<_para_num;i++)
		_res[i]=_best_g[i];

	return(_best_g[_para_num]);
}


