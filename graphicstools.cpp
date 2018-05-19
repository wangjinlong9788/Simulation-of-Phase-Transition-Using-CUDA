//---------------------------------------------------------------------------
// graphicstools.cpp
//---------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <string>

#pragma hdrstop

#include "fileutils.h"
#include "stringutils.h"
#include "graphicstools.h"

//---------------------------------------------------------------------------
#define ABS(a) ((a)<0.0?(-(a)):(a))
//---------------------------------------------------------------------------


colorfunction::colorfunction() 
{
    int i;
    gradidx=0;
    gradID=0;
    grads=new gradient[MAXGRADS];
    for (i=0; i<MAXGRADS; i++) {
        grads[i].ID=i+1;
        grads[i].N=0;
        grads[i].type=1;//RGB
        grads[i].xcols=NULL;
    }
    //set the built-in gradients
    Ngrads=10;
    grads[0].N=17;
    grads[0].xcols=CGrainbow;
    grads[1].N=7;
    grads[1].xcols=CGrainbow2;
    grads[2].N=7;
    grads[2].xcols=CGrainbow3;
    grads[3].N=11;
    grads[3].xcols=CGdarkrainbow;
    grads[4].N=12;
    grads[4].xcols=CGtemperature;
    grads[5].N=13;
    grads[5].xcols=CGtemperature2;
    grads[6].N=12;
    grads[6].xcols=CGthermo;
    grads[7].N=5;
    grads[7].xcols=CGsolar;
    grads[8].N=7;
    grads[8].xcols=CGsunset;
    grads[9].N=10;
    grads[9].xcols=CGneon;
};

colorfunction::~colorfunction() 
{
    int i;
    for(i=0;i<Ngrads;i++)
        if(grads[i].ID&GRAD_USER_LOADED) delete[] grads[i].xcols;
    delete[] grads;
};

//---------------------------------------------------------------------------

//get a RGBA value for 0<=val<=1, using current gradient 
dcolor colorfunction::getgradientcolor(double val,bool periodic)
{
    int n,Ncol;
    dcolor res;
    double dx,dc,f;
    Xcolor *xcol=grads[gradidx].xcols;
    Ncol=grads[gradidx].N;
    
    if(periodic)
    {
        if(val>1.0) val=val-((int) val);
        else if(val<0.0) val=val+1-((int) val);
    }
    
    n=0;
    if(val<=xcol[0].X1) {res.R=xcol[0].R;res.G=xcol[0].G;res.B=xcol[0].B;res.A=xcol[0].A;return res;}
    else if(val>=xcol[Ncol-1].X1) {res.R=xcol[Ncol-1].R;res.G=xcol[Ncol-1].G;res.B=xcol[Ncol-1].B;res.A=xcol[Ncol-1].A;return res;}
 
    
    n=Ncol-1;
    while(val<xcol[n].X1) n--;
    val-=xcol[n].X1;
    dx=val/(xcol[n+1].X1-xcol[n].X1);
 
    f=gnufunc(dx,xcol[n].X2); //standard in linear (3 or negative)
    
    res.R=xcol[n].R;
    dc=xcol[n+1].R-res.R;
    res.R+=f*dc;
 
    res.G=xcol[n].G;
    dc=xcol[n+1].G-res.G;
    res.G+=f*dc;
 
    res.B=xcol[n].B;
    dc=xcol[n+1].B-res.B;
    res.B+=f*dc;
 
    res.A=xcol[n].A;
    dc=xcol[n+1].A-res.A;
    res.A+=f*dc;
 
    cvalidate(res);
    return res; 
};

unsigned int colorfunction::loadgradient(string fn) //returns the ID, does not select
//file format: not defined yet, XML with {R,G,B,A,X1,X2} info (with default A=0,X2=3 and X1=0...1 equidistant.)
{
    return 0; //no implemented
};


unsigned int colorfunction::definegradient(int N,Xcolor *cols) //sets a defined gradient (array cols is used and not copied -> deallocation by user), returns ID, no select
{   
    unsigned int res=0;
    if(Ngrads<MAXGRADS)
    {
        res=GRAD_USER_SET+Ngrads+1;
        grads[Ngrads].N=N;
        grads[Ngrads].ID=res;
        grads[Ngrads].xcols=cols;
        Ngrads++;
    }
    return res;
};


int colorfunction::selectgradient(unsigned int ID) //built-in or user defined, return=0 success
{
    int i;
    for(i=0;i<Ngrads;i++) 
    {
        if(grads[i].ID==ID) {gradidx=i;gradID=ID;return 0;}
    }
    return -1; //invalid ID
};

int colorfunction::selectgradientidx(unsigned int idx) //built-in or user defined, return=0 success
{
    if(idx>=Ngrads) return -1;//invalid idx
    gradidx=idx;
    gradID=grads[idx].ID;
    return 0; 
};

//---------------------------------------------------------------------------

//these are the 37 gnuplot rgbformulae
double inline colorfunction::gnufunc(double val, int ID)
{
  double y;
  if(val<0) val=0;
  else if(val>1.0) val=1.0;

  if(ID<0) {val=1-val;ID=-ID;};
  
  switch(ID)
  {
   case 0:y=0;break;
   case 1:y=0.5;break;
   case 2:y=1;break;
   //case 3:break; //identity
   case 6:y=val*val;
   case 4:y=val*val;break;
   case 5:y=val*val*val;break;
   case 8:y=sqrt(val);
   case 7:y=sqrt(val);break;
   case 9:y=sin(1.570796326794897*val);break;
   case 10:y=cos(1.570796326794897*val);break;
   case 11:y=ABS(val-0.5);break;
   case 12:y=val+val-1;y=y*y;break;
   case 13:y=sin(3.141592653589793*val);break;
   case 14:y=ABS(cos(3.141592653589793*val));break;
   case 15:y=sin(6.283185307179586*val);break;
   case 16:y=cos(6.283185307179586*val);break;
   case 17:y=ABS(sin(6.283185307179586*val));break;
   case 18:y=ABS(cos(6.283185307179586*val));break;
   case 19:y=ABS(sin(12.566370614359172*val));break;
   case 20:y=ABS(cos(12.566370614359172*val));break;
   case 21:y=3*val;break;
   case 22:y=3*val-1;break;
   case 23:y=3*val-2;break;
   case 24:y=ABS(3*val-1);break;
   case 25:y=ABS(3*val-2);break;
   case 26:y=1.5*val-0.5;break;
   case 27:y=1.5*val-1;break;
   case 28:y=ABS(1.5*val-0.5);break;
   case 29:y=ABS(1.5*val-1);break;
   case 30:y=3.125*val-0.78125;break;
   case 31:y=val+val-0.84;break;
   case 32:y=4*val;
           if(y>1) y=-2*val+1.84;
           if(y<0) y=12.5*val-11.5;
           break; 
   case 33:y=ABS(val+val-0.5);break;
   case 34:y=val+val;break;
   case 35:y=val+val-0.5;break;
   case 36:y=val+val-1.0;break;
   default: y=val;
  }
  if(y<0) y=0;
  else if(y>1.0) y=1.0;

  return y;
};
//---------------------------------------------------------------------------

dcolor colorfunction::getgnugradientcolor(double val,int rf,int gf,int bf,int af)
{
    dcolor res;
    if(gf<0) gf=rf;
    if(bf<0) bf=rf;
    res.R=gnufunc(val,rf);
    res.G=gnufunc(val,gf);  
    res.B=gnufunc(val,bf);
    res.A=gnufunc(val,af);
    return res;
};
//---------------------------------------------------------------------------

dcolor colorfunction::RGBtoHSV(dcolor RGB)
{
 dcolor HSV;
 double min, max, delta;
 int m=1;
 HSV.A=RGB.A;
 
 min=max=RGB.R;
 if(RGB.G>max) {max=RGB.G;m=2;} else min=RGB.G;
 if(RGB.B>max) {max=RGB.B;m=3;}
 if(RGB.B<min) min=RGB.B;

 HSV.V = max;				// v
 delta = max - min;
	if( max != 0 )
		HSV.S = delta / max;		// s
	else {
		// r = g = b = 0		// s = 0, v is undefined
		HSV.S = 0;
		HSV.H = 0; //return black
		return HSV;
	}
	if( m==1 )
		HSV.H = ( RGB.G - RGB.B ) / delta;		// between yellow & magenta
	else if( m==2 )
		HSV.H = 2 + ( RGB.B - RGB.R ) / delta;	// between cyan & yellow
	else
		HSV.H = 4 + ( RGB.R - RGB.G ) / delta;	// between magenta & cyan
	HSV.H *= 0.1666666666666666667;				// normalize to 1
    if(HSV.H<0) HSV.H+=1.0;
    
    cvalidate(HSV);
    return HSV;
};
//---------------------------------------------------------------------------

dcolor colorfunction::HSVtoRGB(dcolor HSV)
{
 int i;
 double f, p, q, t;
 dcolor RGB;
 RGB.A=HSV.A;
 

	if( HSV.S == 0 ) {
		// achromatic (grey)
		RGB.R = RGB.G = RGB.B = HSV.V;
		return RGB;
	}		
	i = (int) (6.0*HSV.H);// sector 0 to 5
	f = HSV.H - i;			// factorial part of h
	p = HSV.V * ( 1 - HSV.S );
	q = HSV.V* ( 1 - HSV.S * f );
	t = HSV.V * ( 1 - HSV.S * ( 1 - f ) );
	switch( i ) {
		case 0:
			RGB.R = HSV.V;
			RGB.G = t;
			RGB.B = p;
			break;
		case 1:
			RGB.R = q;
			RGB.G = HSV.V;
			RGB.B = p;
			break;
		case 2:
			RGB.R = p;
			RGB.G = HSV.V;
			RGB.B = t;
			break;
		case 3:
			RGB.R = p;
			RGB.G = q;
			RGB.B = HSV.V;
			break;
		case 4:
			RGB.R = t;
			RGB.G = p;
			RGB.B = HSV.V;
			break;
		default:		// case 5:
			RGB.R = HSV.V;
			RGB.G = p;
			RGB.B = q;
			break;
	}
 
 
 return RGB;
};

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


PXMfile::PXMfile(string fn, int type)
{
    name=fn;
    Lx=1;
    Ly=1;
    Lz=1;
    bits=8;
    maxgray=255;
    BWthres=0.5;
    channels=3;
    ctype=1; //RGB
    frate=25.0;
    fcount=0;
    fop=false;
    settype(type);
    colf=new colorfunction();
};


PXMfile::~PXMfile()
{
    if(fop) FileClose(f);
    delete colf;
};

int PXMfile::writefile(dcolor *col,int N)
{
    string head=composehead();
    fHandle pf;
    unsigned char *buf;
    unsigned int rgba;
    int i,n;
    
    if(N<Lx*Ly*Lz) return -2; //not enough data
    
    N=Lx*Ly*Lz; //we ignore the rest
    
    pf=FileCreate(name+getext());
    if(pf!=-1)
	{
        FileWrite(pf,head.c_str(),head.length());
        switch(ptype)
        {
            case PXM_P0:
                break;
            case PXM_P1:case PXM_P4:case PXM_P2:case PXM_P5:
                buf=new unsigned char[N];
                for(i=0;i<N;i++) buf[i]=(unsigned char) colf->getgrayscale(col[i],0xFF);
                if(ptype==PXM_P5) FileWrite(pf,buf,N);
                else if(ptype==PXM_P4) {n=convertbitarray(buf,N,buf);FileWrite(pf,buf,n);}
                else if(ptype==PXM_P2) writeascii(pf,buf,N);
                else //BW, ascii
                {
                    rgba=(unsigned int) (BWthres*255);
                    for(i=0;i<N;i++) {if(buf[i]>rgba) buf[i]=1;else buf[i]=0;}
                    writeascii(pf,buf,N);
                }
                delete[] buf;
                break;
            case PXM_P3:case PXM_P6:
                buf=new unsigned char[3*N];
                for(i=0;i<N;i++) {
                    rgba=colf->get32bitcolor(col[i]);
                    buf[3*i]=(rgba>>16) & 0xFF;
                    buf[3*i+1]=(rgba>>8) & 0xFF;
                    buf[3*i+2]=rgba & 0xFF;
                }
                if(ptype==PXM_P6) FileWrite(pf,buf,3*N);
                else writeascii(pf,buf,3*N);
                delete[] buf;
                break;
        }
        FileClose(pf);
    }
    else return -1;
    return 0;
};

int PXMfile::writefile(unsigned int *rgba,int N,bool isrgb)
{
    string head=composehead();
    fHandle pf;
    unsigned char *buf;
    unsigned int c;
    int i,j,n;
    
    if(N<Lx*Ly*Lz) return -2; //not enough data
    
    N=Lx*Ly*Lz; //we ignore the rest
    
    pf=FileCreate(name+getext());
    if(pf!=-1)
	{
        FileWrite(pf,head.c_str(),head.length());
        switch(ptype)
        {
            case PXM_P0:
                break;
            case PXM_P1:case PXM_P4:case PXM_P2:case PXM_P5:
                buf=new unsigned char[N];
                for(i=0;i<N;i++) buf[i]=(unsigned char) colf->getgrayscale(rgba[i],0xFF);
                if(ptype==PXM_P5) FileWrite(pf,buf,N);
                else if(ptype==PXM_P4) {n=convertbitarray(buf,N,buf);FileWrite(pf,buf,n);}
                else if(ptype==PXM_P2) writeascii(pf,buf,N);
                else
                {
                    c=(unsigned int) (BWthres*255);
                    for(i=0;i<N;i++) {if(buf[i]>c) buf[i]=1;else buf[i]=0;}
                    writeascii(pf,buf,N);
                }
                delete[] buf;
                break;
            case PXM_P3:case PXM_P6:
                buf=new unsigned char[3*N];
                for(i=0;i<N;i++){
                    c=rgba[i];
                    buf[3*i]=(c>>16) & 0xFF;
                    buf[3*i+1]=(c>>8) & 0xFF;
                    buf[3*i+2]=c & 0xFF;
                }
                if(ptype==PXM_P6) FileWrite(pf,buf,3*N);
                else writeascii(pf,buf,3*N);
                delete[] buf;
                break;
                
        }
        FileClose(pf);
    }
    else return -1;
    return 0;
};

int PXMfile::writefile(double *gray,int N)
{
    string head=composehead();
    fHandle pf;
    unsigned char *buf;
    unsigned int c;
    dcolor dcol;
    int i,j,n;
    
    if(N<Lx*Ly*Lz) return -2; //not enough data
    
    N=Lx*Ly*Lz; //we ignore the rest
    
    
    pf=FileCreate(name+getext());
    if(pf!=-1)
	{
        FileWrite(pf,head.c_str(),head.length());
        switch(ptype)
        {
            case PXM_P0:
                break;
            case PXM_P1:case PXM_P4:case PXM_P2:case PXM_P5:
                buf=new unsigned char[N];
                for(i=0;i<N;i++) {c=(unsigned int) (256*gray[i]);if(c>255) c=255;buf[i]=(c&0xFF);}
                if(ptype==PXM_P5) FileWrite(pf,buf,N);
                else if(ptype==PXM_P4) {n=convertbitarray(buf,N,buf);FileWrite(pf,buf,n);}
                else if(ptype==PXM_P2) writeascii(pf,buf,N);
                else
                {
                    c=(unsigned int) (BWthres*255);
                    for(i=0;i<N;i++) {if(buf[i]>c) buf[i]=1;else buf[i]=0;}
                    writeascii(pf,buf,N);
                }
                delete[] buf;
                break;
            case PXM_P3:case PXM_P6: //use color gradient of colf to make it color
                buf=new unsigned char[3*N];
                for(i=0;i<N;i++){
                    dcol=colf->getgradientcolor(0.0039216*gray[i]);
                    c=colf->get32bitcolor(dcol);
                    buf[3*i]=(c>>16) & 0xFF;
                    buf[3*i+1]=(c>>8) & 0xFF;
                    buf[3*i+2]=c & 0xFF;
                }
                if(ptype==PXM_P6) FileWrite(pf,buf,3*N);
                else writeascii(pf,buf,3*N);
                delete[] buf;
                break;
        }
        FileClose(pf);
    }
    else return -1;
    return 0;
};


int PXMfile::writefile(unsigned char *gray,int N)
{
    string head=composehead();
    fHandle pf;
    unsigned char *buf;
    unsigned int c;
    dcolor dcol;
    int i,j,n;
    
    if(N<Lx*Ly*Lz) return -2; //not enough data
    
    N=Lx*Ly*Lz; //we ignore the rest
    
    
    pf=FileCreate(name+getext());
    if(pf!=-1)
	{
        FileWrite(pf,head.c_str(),head.length());
        switch(ptype)
        {
            case PXM_P0:
                break;
            case PXM_P1:case PXM_P4:
                buf=new unsigned char[N];
                if(ptype==PXM_P4) {n=convertbitarray(gray,N,buf);FileWrite(pf,buf,n);}
                else
                {
                    c=(unsigned int) (BWthres*maxgray);
                    for(i=0;i<N;i++) {if(gray[i]>c) buf[i]=1;else buf[i]=0;}
                    writeascii(pf,buf,N);
                }
                delete[] buf;
                break;
            case PXM_P2:case PXM_P5:
                if(ptype==PXM_P5) FileWrite(pf,gray,N); //fasted method of all
                else writeascii(pf,gray,N);
                break;
            case PXM_P3:case PXM_P6: //use color gradient of colf to make it color
                buf=new unsigned char[3*N];
                for(i=0;i<N;i++){
                    dcol=colf->getgradientcolor(0.0039216*gray[i]);
                    c=colf->get32bitcolor(dcol);
                    buf[3*i]=(c>>16) & 0xFF;
                    buf[3*i+1]=(c>>8) & 0xFF;
                    buf[3*i+2]=c & 0xFF;
                }
                if(ptype==PXM_P6) FileWrite(pf,buf,3*N);
                else writeascii(pf,buf,3*N);
                delete[] buf;
                break;
                
        }
        FileClose(pf);
    }
    else return -1;
    return 0;
};

int PXMfile::startanifile()
{
};

int PXMfile::writeframe(dcolor *col,int N)
{
};

//very slow, binary preferred
void PXMfile::writeascii(fHandle af,unsigned char* buf,int N)
{
    string s;
    int j,i;
    for(j=0;j<Ly*Lz;j++)
    {
        s="";
        for(i=0;i<Lx;i++) {s=s+" "+IntToStr(buf[j*Lx+i]);}
        s=s+"\n";
        FileWrite(af,s.c_str(),s.length());
    }
};

//---------------------------------------------------------------------------
int PXMfile::checkendian()
{
    unsigned int *nuxi;
    char NUXI[]="UNIX"; //0x554E4958
    int endian=-1;
    nuxi = (unsigned int *) NUXI;
    if(*nuxi==0x554E4958) endian=1; //big
    else if(*nuxi==0x58494E55) endian=0; //little
    else if (*nuxi==0x4E555849) endian=2; //middle
    return endian;
}

string PXMfile::composehead()
{string res;
    switch (ptype) {
        case PXM_P1:res="P1\n"+IntToStr(Lx)+ " "+IntToStr(Ly)+"\n";
            break;
        case PXM_P2:res="P2\n"+IntToStr(Lx)+ " "+IntToStr(Ly)+"\n"+IntToStr(maxgray)+"\n";
            break;
        case PXM_P3:res="P3\n"+IntToStr(Lx)+ " "+IntToStr(Ly)+"\n"+IntToStr(maxgray)+"\n";
            break;
        case PXM_P4:res="P4\n"+IntToStr(Lx)+ " "+IntToStr(Ly)+"\n";
            break;
        case PXM_P5: //limited to 8-bit !!!
            if(maxgray>255) maxgray=255;
            res="P5\n"+IntToStr(Lx)+ " "+IntToStr(Ly)+"\n"+IntToStr(maxgray)+"\n";
            break;
        case PXM_P6://limited to 8-bit !!!
            if(maxgray>255) maxgray=255;
            res="P6\n"+IntToStr(Lx)+ " "+IntToStr(Ly)+"\n"+IntToStr(maxgray)+"\n";
            break;    
        default: //P0 format
            res="P0\nSIZE "+IntToStr(Lx)+ " "+IntToStr(Ly);
            if(Lz>1) res=res+IntToStr(Lz);
            res=res+"\nCH "+IntToStr(channels)+"\nBITS "+IntToStr(bits)+"\n";
            if ((channels==1) && (bits>1)) res=res+"MAX " + IntToStr(maxgray)+"\n"; //no check if maxcol<2^bits
            if (binary) res=res+"DATA B\n";else res=res+"DATA A\n";
            if(bits>8)
            {
                if(checkendian()==0) res=res+"ENDIAN L\n"; else res=res+"ENDIAN B\n"; //all non-L or B need to be converted to B
            }
            res=res+"MODEL "+string(ctypes[cmodel])+"\n";
            if(fcount>1) //multi-image or movie file
            {
                res=res+"FRAMES "+IntToStr(fcount)+"\n";
                res=res+"FRAMERATE "+DoubleToStr(frate)+"\n";
            }
            
            res=res+"ENDHDR\n"; //data comes next
            break;
    }
    return res;
};

//output can be equal to input
int PXMfile::convertbitarray(unsigned char *input,int N,unsigned char *output)
{
    int i,j,m;
    int L;
    unsigned char b,t;
    L=(N+7)/8;
    m=(int) (256*BWthres); if(m>255) m=255;
    t=m&0xFF;
    for(i=0;i<L;i++)
    {
        b=0;
        m=(i<<3);
        for (j=0;j<8; j++) {if((m+j<N) && (input[m+j]>t)) b=(b|(1<<j));}
        output[i]=b;
    }
    return L;
};
  
//---------------------------------------------------------------------------
