#include "complexify.h"
#include <math.h>

#include "sundialstypes.h"

#include "nvector_serial.h"


typedef struct {realtype * p;realtype * * P[15][15];realtype * * Jbd[15][15];integertype * pivot[15][15];realtype q4;realtype om;realtype dx;realtype dz;realtype hdco;realtype haco;realtype vdco;} UserDataStruct;
typedef UserDataStruct *UserData;
extern void f_Cplx_(cplx t,N_Vector y_Cplx_Re_,N_Vector y_Cplx_Im_,N_Vector ydot_Cplx_Re_,N_Vector ydot_Cplx_Im_,UserData f_data_Cplx_Re_,UserData f_data_Cplx_Im_);

extern "C"{
void f_C_Cplx(realtype t_Cplx_Re_,realtype t_Cplx_Im_,N_Vector y_Cplx_Re_,N_Vector y_Cplx_Im_,N_Vector ydot_Cplx_Re_,N_Vector ydot_Cplx_Im_,UserData f_data_Cplx_Re_,UserData f_data_Cplx_Im_);
}

extern void f_Cplx_(cplx t,N_Vector y_Cplx_Re_,N_Vector y_Cplx_Im_,N_Vector ydot_Cplx_Re_,N_Vector ydot_Cplx_Im_,UserData f_data_Cplx_Re_,UserData f_data_Cplx_Im_);
void f_C_Cplx(realtype t_Cplx_Re_,realtype t_Cplx_Im_,N_Vector y_Cplx_Re_,N_Vector y_Cplx_Im_,N_Vector ydot_Cplx_Re_,N_Vector ydot_Cplx_Im_,UserData f_data_Cplx_Re_,UserData f_data_Cplx_Im_)
{
	
	f_Cplx_(cplx(t_Cplx_Re_,t_Cplx_Im_),y_Cplx_Re_,y_Cplx_Im_,ydot_Cplx_Re_,ydot_Cplx_Im_,f_data_Cplx_Re_,f_data_Cplx_Im_);

}
void f_Cplx_(cplx t,N_Vector y_Cplx_Re_,N_Vector y_Cplx_Im_,N_Vector ydot_Cplx_Re_,N_Vector ydot_Cplx_Im_,UserData f_data_Cplx_Re_,UserData f_data_Cplx_Im_)
{
	cplx tmp;

	cplx q3;

	cplx c1;

	cplx c2;

	cplx c1dn;

	cplx c2dn;

	cplx c1up;

	cplx c2up;

	cplx c1lt;

	cplx c2lt;

	cplx c1rt;

	cplx c2rt;

	cplx czdn;

	cplx czup;

	cplx hord1;

	cplx hord2;

	cplx horad1;

	cplx horad2;

	cplx qq1;

	cplx qq2;

	cplx qq3;

	cplx qq4;

	cplx rkin1;

	cplx rkin2;

	cplx s;

	cplx vertd1;

	cplx vertd2;

	cplx zdn;

	cplx zup;

	cplx q4coef;

	cplx delz;

	cplx verdco;

	cplx hordco;

	cplx horaco;

	realtype *ydata_Cplx_Re_;
	realtype *ydata_Cplx_Im_;

	realtype *dydata_Cplx_Re_;
	realtype *dydata_Cplx_Im_;

	int jx;

	int jz;

	int idn;

	int iup;

	int ileft;

	int iright;

	UserData data_Cplx_Re_;
	UserData data_Cplx_Im_;

	cplx Q1;

	cplx Q2;

	cplx C3;

	cplx A3;

	cplx A4;

	cplx KH;

	cplx VEL;

	cplx KV0;

	data_Cplx_Re_=(UserData)f_data_Cplx_Re_;
	data_Cplx_Im_=(UserData)f_data_Cplx_Im_;

	ydata_Cplx_Re_=(*((N_VectorContent_Serial)y_Cplx_Re_->content)).data;
	ydata_Cplx_Im_=(*((N_VectorContent_Serial)y_Cplx_Im_->content)).data;

	dydata_Cplx_Re_=(*((N_VectorContent_Serial)ydot_Cplx_Re_->content)).data;
	dydata_Cplx_Im_=(*((N_VectorContent_Serial)ydot_Cplx_Im_->content)).data;

	Q1=cplx(data_Cplx_Re_->p[0],data_Cplx_Im_->p[0]);

	Q2=cplx(data_Cplx_Re_->p[1],data_Cplx_Im_->p[1]);

	C3=cplx(data_Cplx_Re_->p[2],data_Cplx_Im_->p[2]);

	A3=cplx(data_Cplx_Re_->p[3],data_Cplx_Im_->p[3]);

	A4=cplx(data_Cplx_Re_->p[4],data_Cplx_Im_->p[4]);

	KH=cplx(data_Cplx_Re_->p[5],data_Cplx_Im_->p[5]);

	VEL=cplx(data_Cplx_Re_->p[6],data_Cplx_Im_->p[6]);

	KV0=cplx(data_Cplx_Re_->p[7],data_Cplx_Im_->p[7]);

	s=sin(cplx(data_Cplx_Re_->om,data_Cplx_Im_->om)*t);

	if(s>0.0)
	{
		q3=exp(-A3/s);

		tmp=exp(-A4/s),data_Cplx_Re_->q4=tmp.real(),data_Cplx_Im_->q4=tmp.imag();

	}
	else
	{
		q3=0.0;

		tmp=0.0,data_Cplx_Re_->q4=tmp.real(),data_Cplx_Im_->q4=tmp.imag();

	}
	q4coef=cplx(data_Cplx_Re_->q4,data_Cplx_Im_->q4);

	delz=cplx(data_Cplx_Re_->dz,data_Cplx_Im_->dz);

	verdco=cplx(data_Cplx_Re_->vdco,data_Cplx_Im_->vdco);

	hordco=cplx(data_Cplx_Re_->hdco,data_Cplx_Im_->hdco);

	horaco=cplx(data_Cplx_Re_->haco,data_Cplx_Im_->haco);

	for(jz=0;jz<15;jz++)
	{
		zdn=30.0+(jz-0.5)*delz;

		zup=zdn+delz;

		czdn=verdco*exp(0.2*zdn);

		czup=verdco*exp(0.2*zup);

		idn=(jz == 0)?1:-1;

		iup=(jz == 14)?-1:1;

		for(jx=0;jx<15;jx++)
		{
			c1=cplx(ydata_Cplx_Re_[(0+jx*2+jz*2*15)],ydata_Cplx_Im_[(0+jx*2+jz*2*15)]);

			c2=cplx(ydata_Cplx_Re_[(1+jx*2+jz*2*15)],ydata_Cplx_Im_[(1+jx*2+jz*2*15)]);

			qq1=Q1*c1*C3;

			qq2=Q2*c1*c2;

			qq3=q3*C3;

			qq4=q4coef*c2;

			rkin1=-qq1-qq2+2.0*qq3+qq4;

			rkin2=qq1-qq2-qq4;

			c1dn=cplx(ydata_Cplx_Re_[(0+jx*2+(jz+idn)*2*15)],ydata_Cplx_Im_[(0+jx*2+(jz+idn)*2*15)]);

			c2dn=cplx(ydata_Cplx_Re_[(1+jx*2+(jz+idn)*2*15)],ydata_Cplx_Im_[(1+jx*2+(jz+idn)*2*15)]);

			c1up=cplx(ydata_Cplx_Re_[(0+jx*2+(jz+iup)*2*15)],ydata_Cplx_Im_[(0+jx*2+(jz+iup)*2*15)]);

			c2up=cplx(ydata_Cplx_Re_[(1+jx*2+(jz+iup)*2*15)],ydata_Cplx_Im_[(1+jx*2+(jz+iup)*2*15)]);

			vertd1=czup*(c1up-c1)-czdn*(c1-c1dn);

			vertd2=czup*(c2up-c2)-czdn*(c2-c2dn);

			ileft=(jx == 0)?1:-1;

			iright=(jx == 14)?-1:1;

			c1lt=cplx(ydata_Cplx_Re_[(0+(jx+ileft)*2+jz*2*15)],ydata_Cplx_Im_[(0+(jx+ileft)*2+jz*2*15)]);

			c2lt=cplx(ydata_Cplx_Re_[(1+(jx+ileft)*2+jz*2*15)],ydata_Cplx_Im_[(1+(jx+ileft)*2+jz*2*15)]);

			c1rt=cplx(ydata_Cplx_Re_[(0+(jx+iright)*2+jz*2*15)],ydata_Cplx_Im_[(0+(jx+iright)*2+jz*2*15)]);

			c2rt=cplx(ydata_Cplx_Re_[(1+(jx+iright)*2+jz*2*15)],ydata_Cplx_Im_[(1+(jx+iright)*2+jz*2*15)]);

			hord1=hordco*(c1rt-2.0*c1+c1lt);

			hord2=hordco*(c2rt-2.0*c2+c2lt);

			horad1=horaco*(c1rt-c1lt);

			horad2=horaco*(c2rt-c2lt);

			tmp=vertd1+hord1+horad1+rkin1,dydata_Cplx_Re_[(0+jx*2+jz*2*15)]=tmp.real(),dydata_Cplx_Im_[(0+jx*2+jz*2*15)]=tmp.imag();

			tmp=vertd2+hord2+horad2+rkin2,dydata_Cplx_Re_[(1+jx*2+jz*2*15)]=tmp.real(),dydata_Cplx_Im_[(1+jx*2+jz*2*15)]=tmp.imag();

		}
	}
}
