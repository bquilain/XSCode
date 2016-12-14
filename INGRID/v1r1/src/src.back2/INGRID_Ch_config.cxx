#ifndef _INGRID_Ch_config_C
#define _INGRID_Ch_config_C

#include<iostream>
#include<sstream>
#include<fstream>
#include"INGRID_Ch_config.hxx"
using namespace std;

bool INGRID_Ch_config::channel_configuration(int *rmm,int *tfb,int *trip,int *trip_ch,int *mod,int *plane,int *ch,bool *tpl,bool *veto){

  for(int i=0;i<300;i++){
    if(ReadoutId[i][0]==*rmm&&ReadoutId[i][1]==*tfb){
      *mod=PhysId[i][0];
      offdiag = false;
      if     (*mod<7){horizontal=true;vertical=false;}
      else if(*mod>=7&&*mod<14){horizontal=false;vertical=true;}
      else if(*mod==14||*mod==15){offdiag=true;}
      if(*mod==0){full_veto=true;}
      else if(*mod==7){full_veto=true;}
      else if(offdiag){full_veto=true;}
      else{full_veto=false;}


      //read physical channel
      if(tpl_flag[i]){//tracking plane
	*tpl=true;*veto=false;

	*plane=PhysId[i][1];
	*ch=(*trip)*16+(*trip_ch);
	if(*ch>47){return false;}
	//replace physical channel(add 090928 by Otani)
	if(*mod==3){
	  if(*ch<24){
	    *ch=*ch+24;
	  }
	  else{
	    *ch=47-*ch;
	  }

	}
	else{
	  if(*ch<24){
	    *ch=23-*ch;
	  }
	}

	return true;
      }//tracking plane end

      if(!tpl_flag[i]){//VETO plane
	*tpl=false;*veto=true;
	if(PhysId[i][1]==11){//*TFB==11
	  if(*trip==0){
	    *ch=*trip_ch;
	    *plane=3;
	    //replace physical channel(add 090928 by Otani)
	    *ch=21-*ch;
	    return true;
	  }
	  if(*trip==1){
	    if(*trip_ch<6){
	      *ch=(*trip_ch)+16;
	      *plane=3;
	      //replace physical channel(add 090928 by Otani)
	      *ch=21-*ch;
	      return true;
	    }
	    if(*trip_ch==6||*trip_ch==7)return false;
	    if(*trip_ch>=8){
	      *ch=*trip_ch-8;
	      *plane=0;
	      if(horizontal&&(!(full_veto))){return false;}
	      return true;
	    }
	  }
	  if(*trip==2){
	    *ch=*trip_ch+8;
	    *plane=0;
	    if(horizontal&&(!(full_veto))){return false;}
	    if(*ch>21)return false;

	    return true;
	  }
	  if(*trip==3){
	    return false;
	  }
	}//*TFB==11
	if(PhysId[i][1]==12){//*TFB==12
	  if(*trip==0){
	    *ch=*trip_ch;
	    *plane=2;
	    //replace physical channel(add 090928 by Otani)
	    *ch=21-*ch;
	    if(vertical&&(!(full_veto))){return false;}
	    return true;
	  }
	  if(*trip==1){
	    if(*trip_ch<6){
	      *ch=*trip_ch+16;
	      *plane=2;
	      //replace physical channel(add 090928 by Otani)
	      *ch=21-*ch;
	      if(vertical&&(!(full_veto))){return false;}
	      return true;
	    }
	    if(*trip_ch==6||*trip_ch==7)return false;
	    if(*trip_ch>=8){
	      *ch=*trip_ch-8;
	      *plane=1;
	      return true;
	    }
	  }
	  if(*trip==2){
	    *ch=*trip_ch+8;
	    *plane=1;
	    if(*ch>21)return false;
	    return true;
	  }
	  if(*trip==3){
	    return false;
	  }
	}//*TFB==12
      }

    }
  }
  return false;  
}

bool INGRID_Ch_config::active_pln(int mod, int view ,int plane){
  if(plane<11) return true;
  if(mod == 0){
    if     ( view == 0 ){
      if( plane == 13 || plane == 14 )
	return true;
      else
	return false;
    }
    else if( view == 1 ){
      if( plane == 11 || plane == 12 )
	return true;
      else
	return false;
    }
  }
  if(0 < mod && mod <= 6){
    if     ( view == 0 ){
      if( plane == 13 || plane == 14 )
	return true;
      else
	return false;
    }
    else if( view == 1 ){
      if( plane == 12 )
	return true;
      else
	return false;
    }
  }
  if(mod == 7){
    if     ( view == 0 ){
      if( plane == 13 || plane == 14 )
	return true;
      else
	return false;
    }
    else if( view == 1 ){
      if( plane == 11 || plane == 12 )
	return true;
      else
	return false;
    }
  }
  if(7 < mod && mod <= 13){
    if     ( view == 0 ){
      if( plane == 14 )
	return true;
      else
	return false;
    }
    else if( view == 1 ){
      if( plane == 11 || plane == 12 )
	return true;
      else
	return false;
    }


  }
  return false;
}

bool INGRID_Ch_config::active_veto(int *mod,int *plane){
  if(*mod==0)return true;
  if(0<*mod&&*mod<=6){
    if(*plane==0)return false;
    return true;
  }
  if(*mod==7)return true;
  if(7<*mod&&*mod<=13){
    if(*plane==2)return false;
    return true;
  }
  if(*mod==14||*mod==15)
    return true;
  return false;
}

bool INGRID_Ch_config::get_global_pln(int *mod,int *plane, bool *tpl, bool *veto, int *global_pln){
  *global_pln=0;
  int tempmod=*mod;
  int temppln=*plane;
  if(*veto&&!this->active_veto(&tempmod, &temppln))return false;
  for(int nummod=0;nummod<16;nummod++){
    for(int pln=0;pln<11;pln++){
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

bool INGRID_Ch_config::get_global_ch(int *mod, int *plane, int *ch, bool *tpl, bool *veto, int *global_ch){
  //if( *mod<16 && *ch > 23 )return false;
  *global_ch=0;
  int tempmod=*mod;
  int temppln=*plane;
  if(*veto&&!this->active_veto(&tempmod, &temppln))return false;

  for(int pln=0;pln<11;pln++){
    for(int numch=0;numch<48;numch++){
      if((*plane)==pln&&(*ch)==numch&&*tpl)return true;
      (*global_ch)++;
    }//nuch
  }//pln
  for(int pln=0;pln<4;pln++){
    if(this->active_veto(&tempmod, &pln)){
      for(int numch=0;numch<22;numch++){
	if((*plane)==pln&&(*ch)==numch&&*veto)return true;
	(*global_ch)++;
      }//numch
    }//active
  }//pln

  return false;
}


bool INGRID_Ch_config::channel_configuration(int *rmm,int *tfb,int *trip,int *trip_ch,int *mod, int *view, int *plane,int *ch){
  bool tpl_flag;
  bool veto_flag;
  if(!this->channel_configuration(rmm,tfb,trip,trip_ch,mod,plane,ch,&tpl_flag,&veto_flag))return false;
  if(veto_flag){
    *plane=(*plane)+11;
    if(*plane==11||*plane==12){
      *view=1;
      return true;
    }
    else if(*plane==13||*plane==14){
      *view=0;
      return true;
    }
    return false;
  }
  else if(tpl_flag){
    if(*ch>=24){
      *view=0;
      *ch=(*ch)-24;
      return true;
    }
    else{
      *view=1;
      return true;
    }
  }
}

bool INGRID_Ch_config::get_global_ch(int *mod, int *view, int *plane,int *ch, int *global_ch){
  *global_ch=0;
  int temppln=*plane-11;
  int tempmod=*mod;
  if((*plane>10)&&!this->active_veto(&tempmod, &temppln))return false;
  for(int pln=0;pln<11;pln++){
    for(int numch=0;numch<24;numch++){
      if((*plane)==pln&&(*ch)==numch&&(*view)==0)return true;
      (*global_ch)++;
    }//numchX
    for(int numch=0;numch<24;numch++){
      if((*plane)==pln&&(*ch)==numch&&(*view)==1)return true;
      (*global_ch)++;
    }//numchY
  }//pln
  for(int pln=11;pln<15;pln++){
    int temp=pln-11;
    int v;
    if(pln==11||pln==12)v=1;
    else v=0;
    if(this->active_veto(&tempmod, &temp)){
      for(int numch=0;numch<22;numch++){
	if(*view==v&&(*plane)==pln&&(*ch)==numch)return true;
	(*global_ch)++;
      }
    }
  }
  return false;

}

bool INGRID_Ch_config::get_global_ch(int mod, int view, int plane,int ch, int *global_ch){
  *global_ch=0;
  int temppln=plane-11;
  int tempmod=mod;

  if((plane>10)&&!this->active_veto(&tempmod, &temppln))return false;
  for(int pln=0;pln<11;pln++){
    for(int numch=0;numch<24;numch++){
      if((plane)==pln&&(ch)==numch&&(view)==0)return true;
      (*global_ch)++;
    }//numchX
    for(int numch=0;numch<24;numch++){
      if((plane)==pln&&(ch)==numch&&(view)==1)return true;
      (*global_ch)++;
    }//numchY
  }//pln
  for(int pln=11;pln<15;pln++){
    int temp=pln-11;
    int v;
    if(pln==11||pln==12)v=1;
    else v=0;
    if(this->active_veto(&tempmod, &temp)){
      for(int numch=0;numch<22;numch++){
	if(view==v&&(plane)==pln&&(ch)==numch)return true;
	(*global_ch)++;
      }
    }
  }
  return false;

}



#endif
