#ifndef _PM_Ch_config_C
#define _PM_Ch_config_C

#include<iostream>
#include<sstream>
#include<fstream>
#include"PM_Ch_config.hxx"
using namespace std;

bool PM_Ch_config::channel_configuration(int *rmm,int *tfb,int *trip,int *trip_ch,int *mod,int *plane,int *ch,bool *tpl,bool *veto){

  for(int i=0;i<300;i++){
    if(ReadoutId[i][0]==*rmm&&ReadoutId[i][1]==*tfb){
      *mod=PhysId[i][0];

      horizontal=true;vertical=false;full_veto=true;

      //read physical channel
      if(tpl_flag[i]){//tracking plane
	*tpl=true;*veto=false;

	*plane=PhysId[i][1];
	*ch=(*trip)*16+(*trip_ch);
	if(*plane==0){
	  if(*ch>47){return false;}
	    }
	else if(*ch>63){return false;}

	return true;
      }//tracking plane end

      if(!tpl_flag[i]){//VETO plane
	*tpl=false;*veto=true;
 
	if(PhysId[i][1]==18){//*TFB==18
	  if(*trip==0 || *trip==1){
	    *ch=(*trip)*16+(*trip_ch);
	    *plane=0;
	    if(*ch>16){return false;}
	    return true;
	  }
	  if(*trip==2 || *trip==3){
	    *ch=(*trip-2)*16+(*trip_ch);
	    *plane=1;
	    if(*ch>16){return false;}
	    return true;
	  }
	}//*TFB==18

	if(PhysId[i][1]==19){//*TFB==19
	    *ch=(*trip)*16+(*trip_ch);
	    *plane=2;
	    if(*ch>16){return false;}
	    return true;
	}//*TFB==19


	if(PhysId[i][1]==20){//*TFB==20
	    *ch=(*trip)*16+(*trip_ch);
	    *plane=3;
	    if(*ch>16){return false;}
	    return true;
	}//*TFB==20


      }//VETO plane

    }
  }
  return false;  
}

bool PM_Ch_config::active_veto(int *mod,int *plane){
  return true;
}

bool PM_Ch_config::get_global_pln(int *mod,int *plane, bool *tpl, bool *veto, int *global_pln){
  *global_pln=0;
  int tempmod=*mod;
  int temppln=*plane;
  if(*veto&&!this->active_veto(&tempmod, &temppln))return false;
  for(int nummod=15;nummod<16;nummod++){
    for(int pln=0;pln<18;pln++){
      if((*mod)==nummod&&(*plane)==pln&&*tpl)return true;
      (*global_pln)++; 
    }//tpl
    for(int pln=0;pln<4;pln++){
      if((*mod)==nummod&&(*plane)==pln&&*veto)return true;
      if(this->active_veto(&nummod, &pln))(*global_pln)++; 
    }//veto
  }//nummod
  return false;
}

bool PM_Ch_config::get_global_ch(int *mod, int *plane, int *ch, bool *tpl, bool *veto, int *global_ch){
  *global_ch=0;
  int tempmod=*mod;
  int temppln=*plane;
  if(*veto&&!this->active_veto(&tempmod, &temppln))return false;



  for(int pln=0;pln<1;pln++){
    for(int numch=0;numch<48;numch++){
      if((*plane)==pln&&(*ch)==numch&&*tpl)return true;
      (*global_ch)++;
    }//nuch
  }//pln


  for(int pln=1;pln<18;pln++){
    for(int numch=0;numch<64;numch++){
      if((*plane)==pln&&(*ch)==numch&&*tpl)return true;
      (*global_ch)++;
    }//nuch
  }//pln
  for(int pln=0;pln<4;pln++){
    if(this->active_veto(&tempmod, &pln)){
      for(int numch=0;numch<17;numch++){
	if((*plane)==pln&&(*ch)==numch&&*veto)return true;
	(*global_ch)++;
      }//numch
    }//active
  }//pln

  return false;
}


bool PM_Ch_config::channel_configuration(int *rmm,int *tfb,int *trip,int *trip_ch,int *mod, int *view, int *plane,int *ch){
  bool tpl_flag;
  bool veto_flag;
  if(!this->channel_configuration(rmm,tfb,trip,trip_ch,mod,plane,ch,&tpl_flag,&veto_flag))return false;
  if(veto_flag){
    *plane=(*plane)+18;
    if(*plane==18||*plane==20){
      *view=0;
      return true;
    }
    else if(*plane==19||*plane==21){
      *view=1;
      return true;
    }
    return false;
  }
  else if(tpl_flag){
    if(*plane==0)
      {
	if(*ch>=24){
	  *view=1;
	  *ch=*ch-24;
	  return true;
	}
	else{
	  *view=0;
	  return true;
	}
      }
    else{
      if(*ch>=32){
	*view=1;
	  *ch=*ch-32;
	return true;
      }
      else{
	*view=0;
	return true;
      }
    }
  }
}

bool PM_Ch_config::get_global_ch(int *mod, int *view, int *plane,int *ch, int *global_ch){
  *global_ch=0;
  int temppln=*plane-18;//k
  int tempmod=*mod;

  if((*plane>17)&&!this->active_veto(&tempmod, &temppln))return false;

  for(int pln=0;pln<1;pln++){
    for(int numch=0;numch<24;numch++){
      if((*plane)==pln&&(*ch)==numch&&(*view)==0)return true;
      (*global_ch)++;
    }//numchX
    for(int numch=0;numch<24;numch++){
      if((*plane)==pln&&(*ch)==numch&&(*view)==1)return true;
      (*global_ch)++;
    }//numchY
  }//pln

  for(int pln=1;pln<18;pln++){
    for(int numch=0;numch<32;numch++){
      if((*plane)==pln&&(*ch)==numch&&(*view)==0)return true;
      (*global_ch)++;
    }//numchX
    for(int numch=0;numch<32;numch++){
      if((*plane)==pln&&(*ch)==numch&&(*view)==1)return true;
      (*global_ch)++;
    }//numchY
  }//pln



  for(int pln=18;pln<22;pln++){//k
    int temp=pln-18;//k
    int v;
    if(pln==18||pln==20)v=0;//k
    else v=1;
    if(this->active_veto(&tempmod, &temp)){
      for(int numch=0;numch<17;numch++){//k
	if(*view==v&&(*plane)==pln&&(*ch)==numch)return true;
	(*global_ch)++;
      }
    }
  }
  return false;

}

bool PM_Ch_config::get_global_ch(int mod, int view, int plane,int ch, int *global_ch){
  *global_ch=0;
  int temppln=plane-18;
  int tempmod=mod;

  if((plane>17)&&!this->active_veto(&tempmod, &temppln))return false;


  for(int pln=0;pln<1;pln++){
    for(int numch=0;numch<24;numch++){
      if( plane==pln && ch==numch && view==0 )return true;
      (*global_ch)++;
    }//numchX
    for(int numch=0;numch<24;numch++){
      if( plane==pln && ch==numch && view==1)return true;
      (*global_ch)++;
    }//numchY
  }//pln

  for(int pln=1;pln<18;pln++){
    for(int numch=0;numch<32;numch++){
      if((plane)==pln&&(ch)==numch&&(view)==0)return true;
      (*global_ch)++;
    }//numchX
    for(int numch=0;numch<32;numch++){
      if((plane)==pln&&(ch)==numch&&(view)==1)return true;
      (*global_ch)++;
    }//numchY
  }//pln



  for(int pln=18;pln<22;pln++){//k
    int temp=pln-18;//k
    int v;
    if(pln==18||pln==20)v=0;//k
    else v=1;
    if(this->active_veto(&tempmod, &temp)){
      for(int numch=0;numch<17;numch++){//k
	if(view==v&&(plane)==pln&&(ch)==numch)return true;
	(*global_ch)++;
      }
    }
  }
  return false;

}



#endif
