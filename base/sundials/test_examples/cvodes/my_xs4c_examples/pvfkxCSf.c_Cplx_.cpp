#include "complexify.h"
#include <math.h>

#include "sundialstypes.h"    /* definitions of realtype, integertype            */

#include "nvector_parallel.h" /* definitions of type N_Vector, macro N_VDATA     */

#include "sundialsmath.h"     /* contains SQR macro                              */

#include "mpi.h"

#include "mpi_cplx.hpp"

typedef struct {realtype * p;realtype q4;realtype om;realtype dx;realtype dy;realtype hdco;realtype haco;realtype vdco;realtype uext[98];integertype my_pe;integertype isubx;integertype isuby;integertype nvmxsub;integertype nvmxsub2;MPI_Comm comm;} UserDataStruct;
typedef UserDataStruct *UserData;
static void BSend_Cplx_(MPI_Comm comm,int my_pe,integertype isubx,integertype isuby,integertype dsizex,integertype dsizey,realtype *udata_Cplx_Re_,realtype *udata_Cplx_Im_);

static void BRecvPost_Cplx_(MPI_Comm comm,MPI_Request * request,int my_pe,integertype isubx,integertype isuby,integertype dsizex,integertype dsizey,realtype *uext_Cplx_Re_,realtype *uext_Cplx_Im_,realtype *buffer_Cplx_Re_,realtype *buffer_Cplx_Im_);

static void BRecvWait_Cplx_(MPI_Request * request,integertype isubx,integertype isuby,integertype dsizex,realtype *uext_Cplx_Re_,realtype *uext_Cplx_Im_,realtype *buffer_Cplx_Re_,realtype *buffer_Cplx_Im_);

static void ucomm_Cplx_(cplx t,N_Vector u_Cplx_Re_,N_Vector u_Cplx_Im_,UserData data_Cplx_Re_,UserData data_Cplx_Im_);

static void fcalc_Cplx_(cplx t,realtype *udata_Cplx_Re_,realtype *udata_Cplx_Im_,realtype *dudata_Cplx_Re_,realtype *dudata_Cplx_Im_,UserData data_Cplx_Re_,UserData data_Cplx_Im_);

extern void f_Cplx_(cplx t,N_Vector u_Cplx_Re_,N_Vector u_Cplx_Im_,N_Vector udot_Cplx_Re_,N_Vector udot_Cplx_Im_,UserData f_data_Cplx_Re_,UserData f_data_Cplx_Im_);

static void BSend_Cplx_(MPI_Comm comm,int my_pe,integertype isubx,integertype isuby,integertype dsizex,integertype dsizey,realtype *udata_Cplx_Re_,realtype *udata_Cplx_Im_)
{
	cplx tmp;

	int i;

	int ly;

	integertype offsetu;

	integertype offsetbuf;

	realtype bufleft_Cplx_Re_[10];
	realtype bufleft_Cplx_Im_[10];

	realtype bufright_Cplx_Re_[10];
	realtype bufright_Cplx_Im_[10];

	if(isuby!=0)
	{
		MPI_Send_Cplx(udata_Cplx_Re_+0,udata_Cplx_Im_+0,dsizex,11,my_pe-2,0,comm);

	}
	if(isuby!=1)
	{
		offsetu=4*dsizex;

		MPI_Send_Cplx(udata_Cplx_Re_+offsetu,udata_Cplx_Im_+offsetu,dsizex,11,my_pe+2,0,comm);

	}
	if(isubx!=0)
	{
		for(ly=0;ly<5;ly++)
		{
			offsetbuf=ly*2;

			offsetu=ly*dsizex;

			for(i=0;i<2;i++)
			{
				bufleft_Cplx_Re_[(offsetbuf+i)]=udata_Cplx_Re_[(offsetu+i)];
				bufleft_Cplx_Im_[(offsetbuf+i)]=udata_Cplx_Im_[(offsetu+i)];

			}
		}
		MPI_Send_Cplx(bufleft_Cplx_Re_+0,bufleft_Cplx_Im_+0,dsizey,11,my_pe-1,0,comm);

	}
	if(isubx!=1)
	{
		for(ly=0;ly<5;ly++)
		{
			offsetbuf=ly*2;

			offsetu=offsetbuf*5+8;

			for(i=0;i<2;i++)
			{
				bufright_Cplx_Re_[(offsetbuf+i)]=udata_Cplx_Re_[(offsetu+i)];
				bufright_Cplx_Im_[(offsetbuf+i)]=udata_Cplx_Im_[(offsetu+i)];

			}
		}
		MPI_Send_Cplx(bufright_Cplx_Re_+0,bufright_Cplx_Im_+0,dsizey,11,my_pe+1,0,comm);

	}
}
static void BRecvPost_Cplx_(MPI_Comm comm,MPI_Request * request,int my_pe,integertype isubx,integertype isuby,integertype dsizex,integertype dsizey,realtype *uext_Cplx_Re_,realtype *uext_Cplx_Im_,realtype *buffer_Cplx_Re_,realtype *buffer_Cplx_Im_)
{
	cplx tmp;

	integertype offsetue;

	realtype *bufleft_Cplx_Re_;
	realtype *bufleft_Cplx_Im_;

	realtype *bufright_Cplx_Re_;
	realtype *bufright_Cplx_Im_;

	bufleft_Cplx_Re_=buffer_Cplx_Re_;
	bufleft_Cplx_Im_=buffer_Cplx_Im_;

	bufright_Cplx_Re_=buffer_Cplx_Re_+10;
	bufright_Cplx_Im_=buffer_Cplx_Im_+10;

	if(isuby!=0)
	{
		MPI_Irecv_Cplx(uext_Cplx_Re_+2,uext_Cplx_Im_+2,dsizex,11,my_pe-2,0,comm,request+0);

	}
	if(isuby!=1)
	{
		offsetue=86;

		MPI_Irecv_Cplx(uext_Cplx_Re_+offsetue,uext_Cplx_Im_+offsetue,dsizex,11,my_pe+2,0,comm,request+1);

	}
	if(isubx!=0)
	{
		MPI_Irecv_Cplx(bufleft_Cplx_Re_+0,bufleft_Cplx_Im_+0,dsizey,11,my_pe-1,0,comm,request+2);

	}
	if(isubx!=1)
	{
		MPI_Irecv_Cplx(bufright_Cplx_Re_+0,bufright_Cplx_Im_+0,dsizey,11,my_pe+1,0,comm,request+3);

	}
}
static void BRecvWait_Cplx_(MPI_Request * request,integertype isubx,integertype isuby,integertype dsizex,realtype *uext_Cplx_Re_,realtype *uext_Cplx_Im_,realtype *buffer_Cplx_Re_,realtype *buffer_Cplx_Im_)
{
	cplx tmp;

	int i;

	int ly;

	integertype dsizex2;

	integertype offsetue;

	integertype offsetbuf;

	realtype *bufleft_Cplx_Re_;
	realtype *bufleft_Cplx_Im_;

	realtype *bufright_Cplx_Re_;
	realtype *bufright_Cplx_Im_;

	MPI_Status status;

	bufleft_Cplx_Re_=buffer_Cplx_Re_;
	bufleft_Cplx_Im_=buffer_Cplx_Im_;

	bufright_Cplx_Re_=buffer_Cplx_Re_+10;
	bufright_Cplx_Im_=buffer_Cplx_Im_+10;

	dsizex2=dsizex+4;

	if(isuby!=0)
	{
		MPI_Wait(request+0,&status);

	}
	if(isuby!=1)
	{
		MPI_Wait(request+1,&status);

	}
	if(isubx!=0)
	{
		MPI_Wait(request+2,&status);

		for(ly=0;ly<5;ly++)
		{
			offsetbuf=ly*2;

			offsetue=(ly+1)*dsizex2;

			for(i=0;i<2;i++)
			{
				uext_Cplx_Re_[(offsetue+i)]=bufleft_Cplx_Re_[(offsetbuf+i)];
				uext_Cplx_Im_[(offsetue+i)]=bufleft_Cplx_Im_[(offsetbuf+i)];

			}
		}
	}
	if(isubx!=1)
	{
		MPI_Wait(request+3,&status);

		for(ly=0;ly<5;ly++)
		{
			offsetbuf=ly*2;

			offsetue=(ly+2)*dsizex2-2;

			for(i=0;i<2;i++)
			{
				uext_Cplx_Re_[(offsetue+i)]=bufright_Cplx_Re_[(offsetbuf+i)];
				uext_Cplx_Im_[(offsetue+i)]=bufright_Cplx_Im_[(offsetbuf+i)];

			}
		}
	}
}
static void ucomm_Cplx_(cplx t,N_Vector u_Cplx_Re_,N_Vector u_Cplx_Im_,UserData data_Cplx_Re_,UserData data_Cplx_Im_)
{
	cplx tmp;

	realtype *udata_Cplx_Re_;
	realtype *udata_Cplx_Im_;

	realtype *uext_Cplx_Re_;
	realtype *uext_Cplx_Im_;

	realtype buffer_Cplx_Re_[20];
	realtype buffer_Cplx_Im_[20];

	MPI_Comm comm;

	int my_pe;

	integertype isubx;

	integertype isuby;

	integertype nvmxsub;

	integertype nvmysub;

	MPI_Request request[4];

	udata_Cplx_Re_=(*((N_VectorContent_Parallel)u_Cplx_Re_->content)).data;
	udata_Cplx_Im_=(*((N_VectorContent_Parallel)u_Cplx_Im_->content)).data;

	comm=data_Cplx_Re_->comm;

	my_pe=data_Cplx_Re_->my_pe;

	isubx=data_Cplx_Re_->isubx;

	isuby=data_Cplx_Re_->isuby;

	nvmxsub=data_Cplx_Re_->nvmxsub;

	nvmysub=10;

	uext_Cplx_Re_=data_Cplx_Re_->uext;
	uext_Cplx_Im_=data_Cplx_Im_->uext;

	BRecvPost_Cplx_(comm,request,my_pe,isubx,isuby,nvmxsub,nvmysub,uext_Cplx_Re_,uext_Cplx_Im_,buffer_Cplx_Re_,buffer_Cplx_Im_);

	BSend_Cplx_(comm,my_pe,isubx,isuby,nvmxsub,nvmysub,udata_Cplx_Re_,udata_Cplx_Im_);

	BRecvWait_Cplx_(request,isubx,isuby,nvmxsub,uext_Cplx_Re_,uext_Cplx_Im_,buffer_Cplx_Re_,buffer_Cplx_Im_);

}
static void fcalc_Cplx_(cplx t,realtype *udata_Cplx_Re_,realtype *udata_Cplx_Im_,realtype *dudata_Cplx_Re_,realtype *dudata_Cplx_Im_,UserData data_Cplx_Re_,UserData data_Cplx_Im_)
{
	cplx tmp;

	realtype *uext_Cplx_Re_;
	realtype *uext_Cplx_Im_;

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

	cplx cydn;

	cplx cyup;

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

	cplx ydn;

	cplx yup;

	cplx q4coef;

	cplx dely;

	cplx verdco;

	cplx hordco;

	cplx horaco;

	int i;

	int lx;

	int ly;

	int jx;

	int jy;

	integertype isubx;

	integertype isuby;

	integertype nvmxsub;

	integertype nvmxsub2;

	integertype offsetu;

	integertype offsetue;

	cplx Q1;

	cplx Q2;

	cplx C3;

	cplx A3;

	cplx A4;

	cplx KH;

	cplx VEL;

	cplx KV0;

	isubx=data_Cplx_Re_->isubx;

	isuby=data_Cplx_Re_->isuby;

	nvmxsub=data_Cplx_Re_->nvmxsub;

	nvmxsub2=data_Cplx_Re_->nvmxsub2;

	uext_Cplx_Re_=data_Cplx_Re_->uext;
	uext_Cplx_Im_=data_Cplx_Im_->uext;

	Q1=cplx(data_Cplx_Re_->p[0],data_Cplx_Im_->p[0]);

	Q2=cplx(data_Cplx_Re_->p[1],data_Cplx_Im_->p[1]);

	C3=cplx(data_Cplx_Re_->p[2],data_Cplx_Im_->p[2]);

	A3=cplx(data_Cplx_Re_->p[3],data_Cplx_Im_->p[3]);

	A4=cplx(data_Cplx_Re_->p[4],data_Cplx_Im_->p[4]);

	KH=cplx(data_Cplx_Re_->p[5],data_Cplx_Im_->p[5]);

	VEL=cplx(data_Cplx_Re_->p[6],data_Cplx_Im_->p[6]);

	KV0=cplx(data_Cplx_Re_->p[7],data_Cplx_Im_->p[7]);

	offsetu=0;

	offsetue=nvmxsub2+2;

	for(ly=0;ly<5;ly++)
	{
		for(i=0;i<nvmxsub;i++)
		{
			uext_Cplx_Re_[(offsetue+i)]=udata_Cplx_Re_[(offsetu+i)];
			uext_Cplx_Im_[(offsetue+i)]=udata_Cplx_Im_[(offsetu+i)];

		}
		offsetu=offsetu+nvmxsub;

		offsetue=offsetue+nvmxsub2;

	}
	if(isuby==0)
	{
		for(i=0;i<nvmxsub;i++)
		{
			uext_Cplx_Re_[(2+i)]=udata_Cplx_Re_[(nvmxsub+i)];
			uext_Cplx_Im_[(2+i)]=udata_Cplx_Im_[(nvmxsub+i)];

		}
	}
	if(isuby==1)
	{
		offsetu=3*nvmxsub;

		offsetue=6*nvmxsub2+2;

		for(i=0;i<nvmxsub;i++)
		{
			uext_Cplx_Re_[(offsetue+i)]=udata_Cplx_Re_[(offsetu+i)];
			uext_Cplx_Im_[(offsetue+i)]=udata_Cplx_Im_[(offsetu+i)];

		}
	}
	if(isubx==0)
	{
		for(ly=0;ly<5;ly++)
		{
			offsetu=ly*nvmxsub+2;

			offsetue=(ly+1)*nvmxsub2;

			for(i=0;i<2;i++)
			{
				uext_Cplx_Re_[(offsetue+i)]=udata_Cplx_Re_[(offsetu+i)];
				uext_Cplx_Im_[(offsetue+i)]=udata_Cplx_Im_[(offsetu+i)];

			}
		}
	}
	if(isubx==1)
	{
		for(ly=0;ly<5;ly++)
		{
			offsetu=(ly+1)*nvmxsub-4;

			offsetue=(ly+2)*nvmxsub2-2;

			for(i=0;i<2;i++)
			{
				uext_Cplx_Re_[(offsetue+i)]=udata_Cplx_Re_[(offsetu+i)];
				uext_Cplx_Im_[(offsetue+i)]=udata_Cplx_Im_[(offsetu+i)];

			}
		}
	}
	dely=cplx(data_Cplx_Re_->dy,data_Cplx_Im_->dy);

	verdco=cplx(data_Cplx_Re_->vdco,data_Cplx_Im_->vdco);

	hordco=cplx(data_Cplx_Re_->hdco,data_Cplx_Im_->hdco);

	horaco=cplx(data_Cplx_Re_->haco,data_Cplx_Im_->haco);

	s=sin(cplx(data_Cplx_Re_->om,data_Cplx_Im_->om)*t);

	if(s>0.0)
	{
		q3=exp(-A3/s);

		q4coef=exp(-A4/s);

	}
	else
	{
		q3=0.0;

		q4coef=0.0;

	}
	tmp=q4coef,data_Cplx_Re_->q4=tmp.real(),data_Cplx_Im_->q4=tmp.imag();

	for(ly=0;ly<5;ly++)
	{
		jy=ly+isuby*5;

		ydn=30.0+(jy-0.5)*dely;

		yup=ydn+dely;

		cydn=verdco*exp(0.2*ydn);

		cyup=verdco*exp(0.2*yup);

		for(lx=0;lx<5;lx++)
		{
			jx=lx+isubx*5;

			offsetue=(lx+1)*2+(ly+1)*nvmxsub2;

			c1=cplx(uext_Cplx_Re_[offsetue],uext_Cplx_Im_[offsetue]);

			c2=cplx(uext_Cplx_Re_[(offsetue+1)],uext_Cplx_Im_[(offsetue+1)]);

			qq1=Q1*c1*C3;

			qq2=Q2*c1*c2;

			qq3=q3*C3;

			qq4=q4coef*c2;

			rkin1=-qq1-qq2+2.0*qq3+qq4;

			rkin2=qq1-qq2-qq4;

			c1dn=cplx(uext_Cplx_Re_[(offsetue-nvmxsub2)],uext_Cplx_Im_[(offsetue-nvmxsub2)]);

			c2dn=cplx(uext_Cplx_Re_[(offsetue-nvmxsub2+1)],uext_Cplx_Im_[(offsetue-nvmxsub2+1)]);

			c1up=cplx(uext_Cplx_Re_[(offsetue+nvmxsub2)],uext_Cplx_Im_[(offsetue+nvmxsub2)]);

			c2up=cplx(uext_Cplx_Re_[(offsetue+nvmxsub2+1)],uext_Cplx_Im_[(offsetue+nvmxsub2+1)]);

			vertd1=cyup*(c1up-c1)-cydn*(c1-c1dn);

			vertd2=cyup*(c2up-c2)-cydn*(c2-c2dn);

			c1lt=cplx(uext_Cplx_Re_[(offsetue-2)],uext_Cplx_Im_[(offsetue-2)]);

			c2lt=cplx(uext_Cplx_Re_[(offsetue-1)],uext_Cplx_Im_[(offsetue-1)]);

			c1rt=cplx(uext_Cplx_Re_[(offsetue+2)],uext_Cplx_Im_[(offsetue+2)]);

			c2rt=cplx(uext_Cplx_Re_[(offsetue+3)],uext_Cplx_Im_[(offsetue+3)]);

			hord1=hordco*(c1rt-2.0*c1+c1lt);

			hord2=hordco*(c2rt-2.0*c2+c2lt);

			horad1=horaco*(c1rt-c1lt);

			horad2=horaco*(c2rt-c2lt);

			offsetu=lx*2+ly*nvmxsub;

			tmp=vertd1+hord1+horad1+rkin1,dudata_Cplx_Re_[offsetu]=tmp.real(),dudata_Cplx_Im_[offsetu]=tmp.imag();

			tmp=vertd2+hord2+horad2+rkin2,dudata_Cplx_Re_[(offsetu+1)]=tmp.real(),dudata_Cplx_Im_[(offsetu+1)]=tmp.imag();

		}
	}
}
extern "C"{
void f_C_Cplx(realtype t_Cplx_Re_,realtype t_Cplx_Im_,N_Vector u_Cplx_Re_,N_Vector u_Cplx_Im_,N_Vector udot_Cplx_Re_,N_Vector udot_Cplx_Im_,UserData f_data_Cplx_Re_,UserData f_data_Cplx_Im_);
}

extern void f_Cplx_(cplx t,N_Vector u_Cplx_Re_,N_Vector u_Cplx_Im_,N_Vector udot_Cplx_Re_,N_Vector udot_Cplx_Im_,UserData f_data_Cplx_Re_,UserData f_data_Cplx_Im_);
void f_C_Cplx(realtype t_Cplx_Re_,realtype t_Cplx_Im_,N_Vector u_Cplx_Re_,N_Vector u_Cplx_Im_,N_Vector udot_Cplx_Re_,N_Vector udot_Cplx_Im_,UserData f_data_Cplx_Re_,UserData f_data_Cplx_Im_)
{
	
	f_Cplx_(cplx(t_Cplx_Re_,t_Cplx_Im_),u_Cplx_Re_,u_Cplx_Im_,udot_Cplx_Re_,udot_Cplx_Im_,f_data_Cplx_Re_,f_data_Cplx_Im_);

}
void f_Cplx_(cplx t,N_Vector u_Cplx_Re_,N_Vector u_Cplx_Im_,N_Vector udot_Cplx_Re_,N_Vector udot_Cplx_Im_,UserData f_data_Cplx_Re_,UserData f_data_Cplx_Im_)
{
	cplx tmp;

	realtype *udata_Cplx_Re_;
	realtype *udata_Cplx_Im_;

	realtype *dudata_Cplx_Re_;
	realtype *dudata_Cplx_Im_;

	UserDataStruct *data_Cplx_Re_;
	UserDataStruct *data_Cplx_Im_;

	udata_Cplx_Re_=(*((N_VectorContent_Parallel)u_Cplx_Re_->content)).data;
	udata_Cplx_Im_=(*((N_VectorContent_Parallel)u_Cplx_Im_->content)).data;

	dudata_Cplx_Re_=(*((N_VectorContent_Parallel)udot_Cplx_Re_->content)).data;
	dudata_Cplx_Im_=(*((N_VectorContent_Parallel)udot_Cplx_Im_->content)).data;

	data_Cplx_Re_=(UserDataStruct *)f_data_Cplx_Re_;
	data_Cplx_Im_=(UserDataStruct *)f_data_Cplx_Im_;

	ucomm_Cplx_(t,u_Cplx_Re_,u_Cplx_Im_,data_Cplx_Re_,data_Cplx_Im_);

	fcalc_Cplx_(t,udata_Cplx_Re_,udata_Cplx_Im_,dudata_Cplx_Re_,dudata_Cplx_Im_,data_Cplx_Re_,data_Cplx_Im_);

}
