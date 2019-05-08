
#include <math.h>
#include "LM_Process.h"

LM_Process::LM_Process()
{

}

LM_Process::~LM_Process()
{

}

//==============================================================
//Ep=B-f(res)
//磁场信号真实值与估计值之差
//==============================================================
void LM_Process::set_Ep(double *Ep, double *res)
{
	double a,b,c,m,n,p,bt,x,y,z,R;
	int i,j,k;

	bt=res[ object_dim * _object_num ];  //bt当成一个参数，放到res数组最后

	for (i=0;i<_sensor_num;i++)
	{
		x=_L[i][0];  //_L 传感器坐标
		y=_L[i][1];
		z=_L[i][2];

		Ep[i*3+0]=_B[i][0];  //_B 传感器信号
		Ep[i*3+1]=_B[i][1];
		Ep[i*3+2]=_B[i][2];

		for (j=0;j<_object_num;j++)
		{
        
    		k=object_dim*j;

    		a=res[0+k];b=res[1+k];c=res[2+k];
        	m=res[3+k];n=res[4+k];p=res[5+k];

    		R=sqrt((x-a)*(x-a)+(y-b)*(y-b)+(z-c)*(z-c));

    		Ep[i*3+0]-=bt*((3*(m*(x-a)+n*(y-b)+p*(z-c))*(x-a))/(R*R*R*R*R)-m/(R*R*R));
    		Ep[i*3+1]-=bt*((3*(m*(x-a)+n*(y-b)+p*(z-c))*(y-b))/(R*R*R*R*R)-n/(R*R*R));
    		Ep[i*3+2]-=bt*((3*(m*(x-a)+n*(y-b)+p*(z-c))*(z-c))/(R*R*R*R*R)-p/(R*R*R));

		}
	}

}


//==============================================================
//J为f(res)的雅可比矩阵
//==============================================================
#define xa (x-a)
#define yb (y-b)
#define zc (z-c)
void LM_Process::set_J()
{
	double a,b,c,m,n,p,bt,x,y,z,R;
	double tmp1;
	double R3,R5,R7;
	int i,j,k;

	bt=_res[ object_dim * _object_num ];

	for (i=0;i<_sensor_num;i++)
	{
		x=_L[i][0];
		y=_L[i][1];
		z=_L[i][2];

		_J[(i*3+0)][object_dim * _object_num]=0;
		_J[(i*3+1)][object_dim * _object_num]=0;
		_J[(i*3+2)][object_dim * _object_num]=0;

		for (j=0;j<_object_num;j++)
		{
        
    		k=object_dim*j;

    		a=_res[0+k];b=_res[1+k];c=_res[2+k];
        	m=_res[3+k];n=_res[4+k];p=_res[5+k];

    		R=sqrt((x-a)*(x-a)+(y-b)*(y-b)+(z-c)*(z-c));
			R3=R*R*R; R5=R3*R*R; R7=R5*R*R;

			tmp1=3*m*(x-a)+3*n*(y-b)+3*p*(z-c);

			//bx对各个参数求导		
    		_J[i*3+0][k+0]=bt*((-3*m*xa-tmp1-3*m*xa)/R5+5*tmp1*xa*xa/R7);
			_J[i*3+0][k+1]=bt*((-3*n*xa     -3*m*yb)/R5+5*tmp1*xa*yb/R7);
			_J[i*3+0][k+2]=bt*((-3*p*xa     -3*m*zc)/R5+5*tmp1*xa*zc/R7);
			_J[i*3+0][k+3]=bt*(3*xa*xa/R5-1/R3);
			_J[i*3+0][k+4]=bt* 3*yb*xa/R5;
			_J[i*3+0][k+5]=bt* 3*zc*xa/R5;
			_J[i*3+0][object_dim * _object_num]+=tmp1*xa/R5-m/R3;

            //by对各个参数求导
			_J[i*3+1][k+0]=bt*((-3*m*yb     -3*n*xa)/R5+5*tmp1*xa*yb/R7);
			_J[i*3+1][k+1]=bt*((-3*n*yb-tmp1-3*n*yb)/R5+5*tmp1*yb*yb/R7);
			_J[i*3+1][k+2]=bt*((-3*p*yb     -3*n*zc)/R5+5*tmp1*zc*yb/R7);
			_J[i*3+1][k+3]=bt* 3*yb*xa/R5;
			_J[i*3+1][k+4]=bt*(3*yb*yb/R5-1/R3);
			_J[i*3+1][k+5]=bt* 3*yb*zc/R5;
			_J[i*3+1][object_dim * _object_num]+=tmp1*yb/R5-n/R3;

			//bz对各个参数求导
			_J[i*3+2][k+0]=bt*((-3*m*zc     -3*p*xa)/R5+5*tmp1*zc*xa/R7);
			_J[i*3+2][k+1]=bt*((-3*n*zc     -3*p*yb)/R5+5*tmp1*zc*yb/R7);
			_J[i*3+2][k+2]=bt*((-3*p*zc-tmp1-3*p*zc)/R5+5*tmp1*zc*zc/R7);
			_J[i*3+2][k+3]=bt* 3*zc*xa/R5;
			_J[i*3+2][k+4]=bt* 3*zc*yb/R5;
			_J[i*3+2][k+5]=bt*(3*zc*zc/R5-1/R3);
			_J[i*3+2][object_dim * _object_num]+=tmp1*zc/R5-p/R3;

		}
	}


}

//==============================================================
//(A+uI)Sp=G
//利用消去法求解线形方程
//有解返回true  
//否则返回false
//==============================================================
bool LM_Process::solve()
{
	int i,j,k,n=_para_num;

	for (i=0;i<_para_num;i++)
	{
		_Sp[i]=_G[i];
		for (j=0;j<_para_num;j++)
		{
			_A2[i][j]=_A[i][j];
		}
		_A2[i][i]+=_u[i];
	}

	for(k=0;k<n;k++)
	{
		double d;
		int l;
		d=_A2[k][k];
		l=k;
		for(i=k+1;i<n;i++)
			while(fabs(_A2[i][k])>fabs(d))
			{
				d=_A2[i][k];
				l=i;
			};
			if(d==0)
			{
				return false;
			}
			if(l!=k)
			{
				double solt;
				for(j=k;j<n;j++)
				{
					solt=_A2[l][j];
					_A2[l][j]=_A2[k][j];
					_A2[k][j]=solt;
				}
				solt=_Sp[k];
				_Sp[k]=_Sp[l];
				_Sp[l]=solt;
			}
			for(j=k+1;j<n;j++)
				_A2[k][j]=_A2[k][j]/_A2[k][k];  //主对角上元素变成1，所以该行都除一个_A2[k][k]
			_Sp[k]=_Sp[k]/_A2[k][k];

			for(i=k+1;i<n;i++)
				for(j=k+1;j<n;j++)
					_A2[i][j]=_A2[i][j]-_A2[i][k]*_A2[k][j];
			for(i=k+1;i<n;i++)
				_Sp[i]=_Sp[i]-_A2[i][k]*_Sp[k];
	}
	for(i=n-2;i>=0;i--)
	{
		double sum;
		sum=0.0;
		for(j=i+1;j<n;j++)
			sum=sum+_A2[i][j]*_Sp[j];
		_Sp[i]=_Sp[i]-sum;
	}
	
	return true;
}
//==============================================================
//A=JT*J   G=JT*Ep
//==============================================================
void LM_Process::set_A()
{
	int i,j,k;

	for (i=0;i<_para_num;i++)
	{
		for (j=0;j<_para_num;j++)
		{	
			_A[i][j]=0;
			for(k=0;k<_sensor_num*3;k++)
			{
				//A[i,j]+=JT[i,k]*J[k,j]
				//=>A[i,j]+=J[k,i]*J[k,j]
				_A[i][j]+=_J[k][i]*_J[k][j];
		
			}
		}
	}

	for (i=0;i<_para_num;i++)
	{
		_G[i]=0;
		for(k=0;k<_sensor_num*3;k++)
		{
			//G[i]+=JT[i,k]*Ep[k]
			//=>G[i]+=J[k,i]*Ep[k]
			_G[i]+=_J[k][i]*_Ep[k];
	
		}
		
	}

}

//==============================================================
//return ||value||
//==============================================================
double LM_Process::model(double *value, int len)
{
	int i;
	double res=0;
	for (i=0;i<len;i++)
	{
		res+=value[i]*value[i];
	}
	return sqrt(res);

}

//==============================================================
//return SpT*(u*Sp+G)
//==============================================================
double LM_Process::adjust()
{
	int i;
	double res=0;
	for (i=0;i<_para_num;i++)
	{
		res+=(_Sp[i]*_u[i]+_G[i])*_Sp[i];
	}
	return res;
}

//===============================================================
//LM算法计算过程
//===============================================================
const double t=1e-3;
const double Err=1e-10;
const int    KMAX=300;
double LM_Process::Compute()
{
	int i;
	int count=0;//记录迭代次数
	double v=2,rou;
	int stop;//停止标志


	set_J();
	set_Ep(_Ep,_res);
	set_A();

	stop=(model(_G,_para_num)<=Err);

	for (i=0;i<_para_num;i++)
	{
		_u[i]=_A[i][i]*t;
	}

	while((!stop) && (count<KMAX))
	{
		
		count++;
		do
		{
			//方程无解出错
			if(!solve())
			{
				stop=2;
				break;
			}

			//更新值太小则结束
			if(model(_Sp,_para_num)<=Err*model(_res,_para_num))
			{
				stop=3;
			}
			else
			{
				for (i=0;i<_para_num;i++)
					_res_new[i]=_res[i]+_Sp[i];

				for (i=0;i<_object_num;i++)
				{
					double one=sqrt(_res_new[3+i*6]*_res_new[3+i*6]+
						            _res_new[4+i*6]*_res_new[4+i*6]+
									_res_new[5+i*6]*_res_new[5+i*6]);

					_res_new[3+i*6]=_res_new[3+i*6]/one;
					_res_new[4+i*6]=_res_new[4+i*6]/one;
					_res_new[5+i*6]=_res_new[5+i*6]/one;
                }

				set_Ep(_Ep2,_res_new);

				rou=pow(model(_Ep,_sensor_num*3),2)-pow(model(_Ep2,_sensor_num*3),2);
				rou/=adjust();
				if(rou>0)
				{
					for (i=0;i<_para_num;i++)
						_res[i]=_res_new[i];

					set_J();
					//set_Ep(_Ep);
					for (i=0;i<_sensor_num*3;i++)
						_Ep[i]=_Ep2[i];
					set_A();

					if (model(_G,_para_num)<Err)
						stop=4;
					if (pow(model(_Ep,_sensor_num*3),2)<=Err)
						stop=5;

					for (i=0;i<_para_num;i++)
						_u[i]*=(1/3.0 > (1-pow(2*rou-1,3)) ? 1/3.0 : (1-pow(2*rou-1,3)));
					v=2.0;
				}
				else
				{
					for (i=0;i<_para_num;i++)
						_u[i]*=v;
					v*=2.0;
				}//end if(rou>0)
				
			}//end if(model(_Sp,_para_num)<=Err*model(_res,_para_num))

		}while((rou<=0)&&(!stop));

	}//end while

	return model(_Ep,_sensor_num*3);

}
//===============================FILE END=================================
