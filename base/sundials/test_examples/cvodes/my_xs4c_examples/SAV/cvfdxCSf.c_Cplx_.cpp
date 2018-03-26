#include "complexify.h"
#include "sundialstypes.h"

#include "nvector_serial.h"


typedef struct __T139228220{realtype p[3];} *UserData;
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

	cplx y1;

	cplx y2;

	cplx y3;

	cplx yd1;

	cplx yd3;

	UserData data_Cplx_Re_;
	UserData data_Cplx_Im_;

	cplx p1;

	cplx p2;

	cplx p3;

	y1=cplx((*((N_VectorContent_Serial)y_Cplx_Re_->content)).data[0],(*((N_VectorContent_Serial)y_Cplx_Im_->content)).data[0]);

	y2=cplx((*((N_VectorContent_Serial)y_Cplx_Re_->content)).data[1],(*((N_VectorContent_Serial)y_Cplx_Im_->content)).data[1]);

	y3=cplx((*((N_VectorContent_Serial)y_Cplx_Re_->content)).data[2],(*((N_VectorContent_Serial)y_Cplx_Im_->content)).data[2]);

	data_Cplx_Re_=(UserData)f_data_Cplx_Re_;
	data_Cplx_Im_=(UserData)f_data_Cplx_Im_;

	p1=cplx(data_Cplx_Re_->p[0],data_Cplx_Im_->p[0]);

	p2=cplx(data_Cplx_Re_->p[1],data_Cplx_Im_->p[1]);

	p3=cplx(data_Cplx_Re_->p[2],data_Cplx_Im_->p[2]);

	yd1=tmp=-p1*y1+p2*y2*y3,(*((N_VectorContent_Serial)ydot_Cplx_Re_->content)).data[0]=tmp.real(),(*((N_VectorContent_Serial)ydot_Cplx_Im_->content)).data[0]=tmp.imag();

	yd3=tmp=p3*y2*y2,(*((N_VectorContent_Serial)ydot_Cplx_Re_->content)).data[2]=tmp.real(),(*((N_VectorContent_Serial)ydot_Cplx_Im_->content)).data[2]=tmp.imag();

	tmp=-yd1-yd3,(*((N_VectorContent_Serial)ydot_Cplx_Re_->content)).data[1]=tmp.real(),(*((N_VectorContent_Serial)ydot_Cplx_Im_->content)).data[1]=tmp.imag();

}
