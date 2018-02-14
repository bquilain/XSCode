#ifndef _INGRID_Dimension_C
#define _INGRID_Dimension_C

#include<iostream>
#include<sstream>
#include<fstream>
#include"INGRID_Dimension.hh"
#include"Const.hh"

using namespace std;

//------------------------------------------------
void Initialize_INGRID_Dimension(){
  char * cINSTALLREPOSITORY = getenv("INSTALLREPOSITORY");

  TFile* f = TFile::Open(Form("%s/MC/wmmc/src/assembly_checklist.root",cINSTALLREPOSITORY),"READ");
  TTree* t = (TTree*) f->Get("tree");
  double xy1,xy3,xy0;
  int view,pln,ch;
  t -> SetBranchAddress("xy1",&xy1);
  t -> SetBranchAddress("xy3",&xy3);
 t -> SetBranchAddress("view",&view);
  t -> SetBranchAddress("pln",&pln);
  t -> SetBranchAddress("ch",&ch);
  for(int i=0;i<t->GetEntries();i++){
    t -> GetEntry(i);
    if(ch==0){
        xy0 = (xy1+xy3)/2./10. ;//offset -- ML 2017/05/05 i remove the +1.
    }
    if(ch<40){
      position_xy[view][pln][ch] = (xy1+xy3)/2./10.; //mm->cm
      position_xy[view][pln][ch] = position_xy[view][pln][ch] - xy0; //mm->cm
    }
    else{
      // we measured the xy position only for plan scintillators
      position_xy[view][pln][ch] = 0;
    }
  }
  t->Delete();
  f->Close();
  f->Delete();

  f = TFile::Open(Form("%s/MC/wmmc/src/position_module_z.root",cINSTALLREPOSITORY),"READ");
  t = (TTree*) f->Get("tree");
  double z;
  double z0;
  t -> SetBranchAddress("z",&z);
  t -> SetBranchAddress("view",&view);
  t -> SetBranchAddress("pln",&pln);
  for(int i=0;i<t->GetEntries();i++){
    t -> GetEntry(i);
    if(i==0){
        z0 = z/10.;//offset
    }
    for(int ch_temp=0;ch_temp<80;ch_temp++){
    	position_z [view][pln][ch_temp] = z/10.; //mm->cm
    	position_z [view][pln][ch_temp] = position_z [view][pln][ch_temp] - z0;
    }
  }
  t->Delete();
  f->Close();
  f->Delete();
  //  std::cout << "call DIMENSION" << std::endl;

  /*  for(int view=0;view<VIEWMAX;view++){
    for(int pln=0;pln<8;pln++){
      for(int ch=0;ch<CHMAX;ch++){
	cout<<view<<" "<<pln<<" "<<ch<<" "<<position_xy[view][pln][ch]<<" "<<position_z[view][pln][ch]<<endl;
      }
    }
    }*/
}


bool INGRID_Dimension::get_pos(int mod, int pln, int ch, bool tpl, bool veto, double *posxy, double *posz){
  if(tpl){
    if(mod!=16){
      *posxy = (ScintiWidth)*ch;
      *posz = (PlnThick+IronThick)*pln;
    }
    else{
      if(pln==0){
	*posxy = (ScintiWidth)*ch;;
	*posz = 0;
      }
      else{
	if(ch<8) *posxy = (ScintiWidth)*ch;
	else if(ch<24)*posxy = (ScintiWidth)*8+(ScibarWidth)*(ch-8);
	else *posxy = (ScintiWidth)*(ch-8);
	*posz = (PlnThick_PM)*(pln-1)+PlnThick_front_PM;
      }
    }
  }
  if(veto){
    /*
    if(pln==11)*posxy =;
    if(pln==12)*posxy =;
    if(pln==13)*posxy =;
    if(pln==14)*posxy =;
    */
    *posxy=0;
    *posz = (ScintiWidth)*ch;
  }
  return true;
}

bool INGRID_Dimension::get_posXY(int mod, int view, int pln, int ch, double *posxy, double *posz){
  // modified ML 2017/06/27 to have position at the center of the bars for PM/INGRID
  if(mod!=16){
    *posxy   = (ScintiWidth)*(ch+.5);//ML
    if( pln <= 10 ){
      if(view == 1) 
	*posz    = ( PlnThick + IronThick ) * pln + ScintiThick/2.;//ML
      else if(view == 0)
	*posz    = ( PlnThick + IronThick ) * pln + ScintiThick + ScintiThick/2.;//ML
      return true;
    }
    else if(pln >= 11){
      this -> get_posVeto( mod, view, pln, ch, posxy, posz );
    }
  }
  else{
    if(pln==0){
      *posxy   = (ScintiWidth)*(ch+.5);//ML
      if(view==0) *posz = ScintiThick/2.;//ML
      else *posz=PlnDist_PM+ScintiThick/2.;//ML
    }
    else if( pln <= 17 ){
      if(ch<8) *posxy = (ScintiWidth)*(ch+0.5);//ML
      else if(ch<24)*posxy = (ScintiWidth)*8+(ScibarWidth)*(ch-8+.5);//ML
      else *posxy = (ScintiWidth)*(ch-8+.5);//ML
      if(view == 0) *posz = (PlnThick_PM)*(pln-1)+PlnThick_front_PM+ScintiThick/2.;//ML
      else *posz = (PlnThick_PM)*(pln-1)+PlnThick_front_PM+PlnDist_PM+ScintiThick/2.;//ML
      return true;
    }
    else if(pln >= 18){
      this -> get_posVeto( mod, view, pln, ch, posxy, posz );
    }
  }
}

bool INGRID_Dimension::get_posVeto(int mod, int view, int pln, int ch, double *posxy, double *posz){
  if(mod!=16){
  if(pln<=10||pln>=15) return false;

  *posz    = VetoStartZ     + ScintiWidth * ch;
  if(pln==11){//############ Right  VETO ################
    //*posz  = VetoOffsetZX + ScintiWidth*ch;
    *posxy = VetoOffsetRight;
  }
  if(pln==12){//############ Left   VETO ################
    //*posz  = VetoOffsetZX + ScintiWidth*ch;
    *posxy = VetoOffsetLeft;
  }
  if(pln==13){//############ Bottom VETO ################
    //*posz  = VetoOffsetZY + ScintiWidth*ch;
    *posxy = VetoOffsetBottom;
  }
  if(pln==14){//############ Up VETO     ################
    //*posz  = VetoOffsetZY + ScintiWidth*ch;
    *posxy = VetoOffsetUp;
  }
  return true;
  }
  else{
  if(pln<=17||pln>=22) return false;
     *posz  = VetoStartZ_PM + ScintiWidth*ch;
    if(pln==21){//Right VETO
      *posxy = VetoOffsetRight_PM;
      return true;
    }
    else if(pln==19){//Left VETO
      *posxy = VetoOffsetLeft_PM;
      return true;
    }
    else if(pln==18){//Bottom VETO
      *posxy = VetoOffsetBottom_PM;
      return true;
    }
    else if(pln==20){//Up VETO
      *posxy = VetoOffsetUp_PM;
      return true;
    }
  }

}

double  Wi = 0.5;

bool INGRID_Dimension::get_expch(int mod, int pln, int *ch, bool tpl, bool veto, double a, double b)
{
  if(mod!=16){
  if(tpl){
    double expz=pln*(PlnThick+IronThick);
    double expxy=expz*a+b;
    for(int numch=0;numch<48;numch++){
      double diffxy=expxy-numch*ScintiWidth;
      if(-Wi*ScintiWidth<=diffxy&&diffxy<Wi*ScintiWidth){
	*ch=numch;
	return true;
      }
    }
    return false;
  }

  if(veto){
    if(pln==0){//Right VETO
      double expz=(VetoOffsetRight-b)/a;
      for(int numch=0;numch<22;numch++){
	double diffxy=expz-numch*ScintiWidth+VetoOffsetZY;
	if(-Wi*ScintiWidth<=diffxy&&diffxy<Wi*ScintiWidth){
	  *ch=numch;
	  return true;
	}
      }
    }
    else if(pln==1){//LEFT VETO
      double expz=(VetoOffsetLeft-b)/a;
      for(int numch=0;numch<22;numch++){
	double diffxy=expz-numch*ScintiWidth+VetoOffsetZY;
	if(-Wi*ScintiWidth<=diffxy&&diffxy<Wi*ScintiWidth){
	  *ch=numch;
	  return true;
	}
      }
    }
    else if(pln==2){//Bottom VETO
      double expz=(VetoOffsetBottom-b)/a;
      for(int numch=0;numch<22;numch++){
	double diffxy=expz-numch*ScintiWidth+VetoOffsetZX;
	if(-Wi*ScintiWidth<=diffxy&&diffxy<Wi*ScintiWidth){
	  *ch=numch;
	  return true;
	}
      }
    }
    else if(pln==3){//UP VETO
      double expz=(VetoOffsetUp-b)/a;
      for(int numch=0;numch<22;numch++){
	double diffxy=expz-numch*ScintiWidth+VetoOffsetZX;
	if(-Wi*ScintiWidth<=diffxy&&diffxy<Wi*ScintiWidth){
	  *ch=numch;
	  return true;
	}
      }
    }

    return false;
  }
  }
  else{

  if(tpl){
    double expz=pln*(PlnThick+IronThick);
    double expxy=expz*a+b;
    for(int numch=0;numch<48;numch++){
      double diffxy=expxy-numch*ScintiWidth;
      if(-0.5*ScintiWidth<=diffxy&&diffxy<0.5*ScintiWidth){
	*ch=numch;
	return true;
      }
    }
    return false;
  }

  if(veto){
    if(pln==3){//Right VETO
      double expz=(VetoOffsetRight_PM-b)/a;
      for(int numch=0;numch<17;numch++){
	double diffxy=expz-numch*ScintiWidth+VetoOffsetZY_PM;
	if(-0.5*ScintiWidth<=diffxy&&diffxy<0.5*ScintiWidth){
	  *ch=numch;
	  return true;
	}
      }
    }
    else if(pln==1){//LEFT VETO
      double expz=(VetoOffsetLeft_PM-b)/a;
      for(int numch=0;numch<17;numch++){
	double diffxy=expz-numch*ScintiWidth+VetoOffsetZY_PM;
	if(-0.5*ScintiWidth<=diffxy&&diffxy<0.5*ScintiWidth){
	  *ch=numch;
	  return true;
	}
      }
    }
    else if(pln==0){//Bottom VETO
      double expz=(VetoOffsetBottom_PM-b)/a;
      for(int numch=0;numch<17;numch++){
	double diffxy=expz-numch*ScintiWidth+VetoOffsetZX_PM;
	if(-0.5*ScintiWidth<=diffxy&&diffxy<0.5*ScintiWidth){
	  *ch=numch;
	  return true;
	}
      }
    }
    else if(pln==2){//UP VETO
      double expz=(VetoOffsetUp_PM-b)/a;
      for(int numch=0;numch<17;numch++){
	double diffxy=expz-numch*ScintiWidth+VetoOffsetZX_PM;
	if(-0.5*ScintiWidth<=diffxy&&diffxy<0.5*ScintiWidth){
	  *ch=numch;
	  return true;
	}
      }
   }

    return false;
  }
  }

}

bool INGRID_Dimension::get_expchXY(int mod, int view, int pln, int *ch, double a, double b){
  int temp=-777;
  if(mod!=16){
  if(pln>=11){//VETO plane
    int veto = pln - 11;
    bool flag = this->get_expch(mod, veto, &temp, 0, 1, a, b);
    *ch = temp;
    return flag;
  }
  else {//Tracking plane
    bool flag = this->get_expch(mod, pln   , &temp, 1, 0, a, b);
    if(temp>23)return false;
    if(temp<0)return false;
    *ch = temp;
    return flag;
  }
  }
  else{
  if(pln>=18){//VETO plane
    int veto = pln - 18;
    bool flag = this->get_expch(mod, veto, &temp, 0, 1, a, b);
    *ch = temp;
    return flag;
  }
  else {//Tracking plane
    bool flag = this->get_expch(mod, pln   , &temp, 1, 0, a, b);
    if(temp>31)return false;
    if(temp<0)return false;
    *ch = temp;
    return flag;
  }
  }
}




//added for prototype of WAGASCI
bool INGRID_Dimension::get_pos_loli(int mod, int view, int pln, int ch, int grid, double *posx, double *posy, double *posz){

	////Based on design
	//double x=0,y=0,z=0;
	//if(view==0){
	//	x = 0;
	//	if(grid==0){
	//		y = loli_offsetxy      + loli_scinti_width/2. + loli_scinti_width * ch; //offset + scinti width + loop
	//		z = loli_firstoffset_z + loli_scinti_thick/2. + loli_gap * pln;		//offset + scinti width + loop
	//	}
	//	else if(grid==1){
	//		y = loli_offsetxy_grid + loli_cutwidth/2.     + loli_cutgap * ch;	//offset + cut width + loop
	//		z = loli_firstoffset_z + loli_scinti_thick    + loli_scinti_width/2. + loli_gap * pln; //offset + offset + scinti width + loop
	//	}
	//	else if(grid==2){
	//		y = loli_offsetxy_grid + loli_cutwidth/2.     + loli_cutgap * ch;	//offset + cut width + loop
	//		z = loli_firstoffset_z + loli_scinti_thick    + loli_offset_hv + loli_scinti_width/2. + loli_gap * pln; //offset + offset + offset + scinti width + loop
	//	}
	//	else return false;
	//}
	//else if(view==1){
	//	y = 0;
	//	if(grid==0){
	//		x = loli_offsetxy      + loli_scinti_width/2. + loli_scinti_width * ch;  //offset + scinti width + loop
	//		z = loli_firstoffset_z + loli_offset_hv + loli_scinti_thick/2. + loli_gap * pln;		//offset + offset + scinti width + loop
	//	}
	//	else if(grid==1){
	//		x = loli_offsetxy_grid + loli_cutwidth/2.     + loli_cutgap * ch;	//offset + cut width + loop
	//		z = loli_firstoffset_z + loli_scinti_thick    +  loli_scinti_width/2. + loli_gap * pln; //offset + offset + offset + scinti width + loop
	//	}
	//	else if(grid==2){
	//		x = loli_offsetxy_grid + loli_cutwidth/2.     + loli_cutgap * ch;	//offset + cut width + loop
	//		z = loli_firstoffset_z + loli_scinti_thick    + loli_offset_hv +  loli_scinti_width/2. + loli_gap * pln; //offset + offset + offset + scinti width + loop
	//	}
	//	else return false;
	//}
	//else return false;
	////std::cout << "INGRID_Dimension " << x << " " << y << " " << z << "\n";
	//*posx = x;
	//*posy = y;
	//*posz = z;
	//return true;


  
	//Based on measurement
	double x=0,y=0,z=0;
	if(view==0){
		x = 0;
		if(grid==0){
			y = loli_offsetxy      + loli_scinti_width/2. + position_xy[view][pln][ch]; //offset + scinti width + loop
			z = loli_firstoffset_z + loli_scinti_thick/2. + position_z [view][pln][ch]; //offset + scinti width + loop
		}
		else if(grid==1){
			y = loli_offsetxy_grid + loli_cutwidth/2.     + loli_cutgap * ch;	//offset + cut width + loop
			z = loli_firstoffset_z + loli_scinti_thick    + loli_scinti_width/2. + position_z [0][pln][ch]; //offset + offset + scinti width + loop
		}
		else if(grid==2){
			y = loli_offsetxy_grid + loli_cutwidth/2.     + loli_cutgap * ch;	//offset + cut width + loop
			z = loli_firstoffset_z + loli_scinti_thick    + loli_scinti_width/2. + position_z [1][pln][ch]; //offset + offset + offset + scinti width + loop
		}
		else return false;
	}
	else if(view==1){
		y = 0;
		if(grid==0){
			x = loli_offsetxy      + loli_scinti_width/2. + position_xy[view][pln][ch];  //offset + scinti width + loop
			z = loli_firstoffset_z + loli_scinti_thick/2. + position_z [view][pln][ch];  //offset + offset + scinti width + loop
		}
		else if(grid==1){
			x = loli_offsetxy_grid + loli_cutwidth/2.     + loli_cutgap * ch;	//offset + cut width + loop
			z = loli_firstoffset_z + loli_scinti_thick    + loli_scinti_width/2. + position_z [0][pln][ch]; //offset + offset + offset + scinti width + loop
		}
		else if(grid==2){
			x = loli_offsetxy_grid + loli_cutwidth/2.     + loli_cutgap * ch;	//offset + cut width + loop
			z = loli_firstoffset_z + loli_scinti_thick    + loli_scinti_width/2. + position_z [1][pln][ch]; //offset + offset + offset + scinti width + loop
		}
		else return false;
	}
	else return false;
	//std::cout << "INGRID_Dimension " << x << " " << y << " " << z << "\n";
	*posx = x;
	*posy = y;
	*posz = z;
	return true;

}


bool INGRID_Dimension::get_grid_loli(int mod, int view, int pln, int ch, int *grid, int *gridch){
	if(ch>=0 && ch<40){
		*grid=0;
		*gridch=ch;
	}
	else if(ch>=40 && ch < 60){
		*grid=1;
		*gridch=ch-40;
	}
	else if(ch>=60 && ch < 80){
		*grid=2;
		*gridch=ch-60;
	}
	else{ 
		std::cout << "INGRID_Dimension  chnum exceed 80 or less than 0" << "\n";
		return false;
	}

	return true;
}

bool INGRID_Dimension::get_pos_loli(int mod, int view, int pln, int ch, double *posx, double *posy, double *posz){
	int grid,gridch;
	double x,y,z;
	this->get_grid_loli(mod, view, pln, ch, &grid, &gridch);
	this->get_pos_loli (mod, view, pln, gridch,  grid, &x, &y, &z);
	*posx=x;
	*posy=y;
	*posz=z;
}


bool INGRID_Dimension::get_pos_loli_xy(int mod, int view, int pln, int ch, double *posxy, double *posz){
	int grid,gridch;
	double x,y,z;
	this->get_grid_loli(mod, view, pln, ch, &grid, &gridch);
	this->get_pos_loli (mod, view, pln, gridch,  grid, &x, &y, &z);
	if(view==0){
		*posxy=y;
	}
	else if(view==1){
		*posxy=x;
	}
	*posz =z;
}



bool INGRID_Dimension::get_loli_gridcell_id(int mod, int view, int pln, int ch, double posx, double posy, double posz,
				int* gridcell_id_x1, int* gridcell_id_x2, int* gridcell_id_y1, int* gridcell_id_y2){	
	int grid, gridch;
	double posx_ch,posy_ch,posz_ch;
	this->get_grid_loli(mod, view, pln, ch, &grid, &gridch);
	if(view==0){
		if(grid==0){
			*gridcell_id_x1= (gridch+1)/2;
			for(int i=0;i<20;i++){
				this->get_pos_loli(mod, 1, pln, i, 1, &posx_ch, &posy_ch, &posz_ch);
				//std::cout << posx_ch << " " << posx << std::endl;
				if(posx_ch > posx){
					*gridcell_id_y1 = i;
					break;
				}
				if(i==19)*gridcell_id_y1 = 20;
			}
		}
		else if(grid==1){
			*gridcell_id_x1= gridch;
			*gridcell_id_x2= gridch+1;
			for(int i=0;i<20;i++){
				this->get_pos_loli(mod, 1, pln, i, 1, &posx_ch, &posy_ch, &posz_ch);
				if(posx_ch > posx){
					*gridcell_id_y1 = i;
					break;
				}
				if(i==19)*gridcell_id_y1 = 20;
			}
		}
		else if(grid==2){
			*gridcell_id_x1= gridch+21;
			*gridcell_id_x2= gridch+21+1;
			for(int i=0;i<20;i++){
				this->get_pos_loli(mod, 1, pln, i, 1, &posx_ch, &posy_ch, &posz_ch);
				if(posx_ch > posx){
					*gridcell_id_y1 = i+21;
					break;
				}
				if(i==19)*gridcell_id_y1 = 41;
			}
		}
	}

	else if(view==1){
		if(grid==0){
			*gridcell_id_y1= (gridch+1)/2 + 21;
			for(int i=0;i<20;i++){
				this->get_pos_loli(mod, 0, pln, i, 1, &posx_ch, &posy_ch, &posz_ch);
				if(posy_ch > posy){
					*gridcell_id_x1 = i+21;
					break;
				}
				if(i==19)*gridcell_id_x1 = 41;
			}
		}
		else if(grid==1){
			*gridcell_id_y1= gridch;
			*gridcell_id_y2= gridch+1;
			for(int i=0;i<20;i++){
				this->get_pos_loli(mod, 0, pln, i, 1, &posx_ch, &posy_ch, &posz_ch);
				if(posy_ch > posy){
					*gridcell_id_x1 = i;
					break;
				}
				if(i==19)*gridcell_id_x1 = 20;
			}
		}
		else if(grid==2){
			*gridcell_id_y1= gridch+21;
			*gridcell_id_y2= gridch+21+1;
			for(int i=0;i<20;i++){
				this->get_pos_loli(mod, 0, pln, i, 1, &posx_ch, &posy_ch, &posz_ch);
				if(posy_ch > posy){
					*gridcell_id_x1 = i+21;
					break;
				}
				if(i==19)*gridcell_id_x1 = 41;
			}
		}
	}

}

#endif
