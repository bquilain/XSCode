#ifndef _INGRID_Dimension_C
#define _INGRID_Dimension_C

#include<iostream>
#include<sstream>
#include<fstream>
#include"INGRID_Dimension.hxx"

using namespace std;

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
        xy0 = (xy1+xy3)/2./10. ;//offset  -- ML 2017/05/05 i remove the +1.
    }
    if(ch<40){
      position_xy[view][pln][ch] = (xy1+xy3)/2./10.; //mm->cm
      position_xy[view][pln][ch] = position_xy[view][pln][ch] - xy0; //mm->cm
    }
    else{
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
    return true;
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
    return false;
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
    return false;
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
	// based on design 
	/*
	if(view==0){
		x = 0;
		if(grid==0){
			y = loli_offsetxy      + loli_scinti_width/2. + loli_scinti_width * ch; //offset + scinti width + loop
			z = loli_firstoffset_z + loli_scinti_thick/2. + loli_gap * pln;		//offset + scinti width + loop
		}
		else if(grid==1){
			y = loli_offsetxy_grid + loli_cutwidth/2.     + loli_cutgap * ch;	//offset + cut width + loop
			z = loli_firstoffset_z + loli_scinti_thick    + loli_scinti_width/2. + loli_gap * pln; //offset + offset + scinti width + loop
		}
		else if(grid==2){
			y = loli_offsetxy_grid + loli_cutwidth/2.     + loli_cutgap * ch;	//offset + cut width + loop
			z = loli_firstoffset_z + loli_scinti_thick    + loli_offset_hv + loli_scinti_width/2. + loli_gap * pln; //offset + offset + offset + scinti width + loop
		}
		else return false;
	}
	else if(view==1){
		y = 0;
		if(grid==0){
			x = loli_offsetxy      + loli_scinti_width/2. + loli_scinti_width * ch;  //offset + scinti width + loop
			z = loli_firstoffset_z + loli_offset_hv + loli_scinti_thick/2. + loli_gap * pln;		//offset + offset + scinti width + loop
		}
		else if(grid==1){
			x = loli_offsetxy_grid + loli_cutwidth/2.     + loli_cutgap * ch;	//offset + cut width + loop
			z = loli_firstoffset_z + loli_scinti_thick    +  loli_scinti_width/2. + loli_gap * pln; //offset + offset + offset + scinti width + loop
		}
		else if(grid==2){
			x = loli_offsetxy_grid + loli_cutwidth/2.     + loli_cutgap * ch;	//offset + cut width + loop
			z = loli_firstoffset_z + loli_scinti_thick    + loli_offset_hv +  loli_scinti_width/2. + loli_gap * pln; //offset + offset + offset + scinti width + loop
		}
		else return false;
	}
	else return false;
	//std::cout << "INGRID_Dimension " << x << " " << y << " " << z << "\n";
	*posx = x;
	*posy = y;
	*posz = z;
	return true;
	*/


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
  if(ch>=0 && ch<40){ // also valid for mod != 15
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
	this->get_grid_loli(mod, view, pln, ch, &grid, &gridch);
	this->get_pos_loli (mod, view, pln, gridch,  grid, posx, posy, posz);
	return true;
}



bool INGRID_Dimension::get_pos_loli_xy(int mod, int view, int pln, int ch, double *posxy, double *posz){
  if(mod==15){
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
  else
    get_posXY(mod,view,pln,ch,posxy,posz);
  
  return true;
}


bool INGRID_Dimension::get_pos_loli_xy(int mod, int view, int pln, int gridch, int grid, double *posxy, double *posz){
	double x,y,z;
	this->get_pos_loli (mod, view, pln, gridch,  grid, &x, &y, &z);
	if(view==0){
		*posxy=y;
	}
	else if(view==1){
		*posxy=x;
	}
	*posz =z;
	return true;
}





////for Lolirecon/////////////

bool INGRID_Dimension::get_reconplnch_loli(int mod, int view, int pln, int ch, int axis, int* reconpln, int* reconch){
	if(mod<=14){
		*reconpln= pln;
		*reconch = ch;
	}
	else if(mod==15){
		int grid,gridch;
		this->get_grid_loli(mod, view, pln, ch, &grid, &gridch);
		if(axis==0){
			if(view==0){
				if(grid==0){
					*reconpln=pln*3;
					*reconch =gridch;
				}
				else if(grid==1){
					*reconpln=pln*3+1;
					*reconch =gridch;
				}
				else if(grid==2){
					*reconpln=pln*3+2;
					*reconch =gridch;
				}
				else return false;
			}
			else if(view==1){
				if(grid==0){
					*reconpln=pln*3+1;
					*reconch =gridch;
				}
				else if(grid==1){
					*reconpln=pln*3;
					*reconch =gridch;
				}
				else if(grid==2){
					*reconpln=pln*3+2;
					*reconch =gridch;
				}
				else return false;
			}
			else return false;
		}
		else if(axis==1){
			if(view==0){
				if(grid==0){
					if(gridch%2==0)*reconpln=(int) (gridch*1.5);
					else if(gridch%2==1)*reconpln=(int) (gridch*1.5+0.5);
					else return false;
					*reconch =pln;
				}
				else if(grid==1){
					*reconpln=gridch*3+1;
					*reconch =pln*2;
				}
				else if(grid==2){
					*reconpln=gridch*3+1;
					*reconch =pln*2+1;
				}
				else return false;
			}
			else if(view==1){
				if(grid==0){
					if(gridch%2==0)*reconpln=(int) (gridch*1.5);
					else if(gridch%2==1)*reconpln=(int) (gridch*1.5+0.5);
					*reconch =pln;
				}
				else if(grid==1){
					*reconpln=gridch*3+1;
					*reconch =pln*2;
				}
				else if(grid==2){
					*reconpln=gridch*3+1;
					*reconch =pln*2+1;
				}
				else return false;
			}
			else return false;
		}
		else return false;
	}
	else if(mod==16){
		*reconpln= pln;
		*reconch = ch;
	}
	else return false;
	
	return true;
}



bool INGRID_Dimension::get_plnch_fromrecon_loli(int mod, int view, int reconpln, int reconch, int axis, int* pln, int* gridch, int* grid){
	if(axis==0){
		*pln=reconpln/3;
		*gridch=reconch;
		if(view==0){
			if(reconpln%3==0){
				*grid=0;
			}
			else if(reconpln%3==1){
				*grid=1;
			}
			else if(reconpln%3==2){
				*grid=2;
			}
		}
		else if(view==1){
			if(reconpln%3==0){
				*grid=1;
			}
			else if(reconpln%3==1){
				*grid=0;
			}
			else if(reconpln%3==2){
				*grid=2;
			}
		}
	}
	else if(axis==1){
		if(reconpln%3==1){
			*pln=reconch/2;
			*gridch=reconpln/3;
			if(reconch%2==0)     *grid=1;
			else if(reconch%2==1)*grid=2;
		}
		else{
			*pln=reconch;
			*gridch = (int) (reconpln/1.5);
			*grid=0;
		}
	}
	
	return true;
}



int INGRID_Dimension::get_chmax_loli(int mod, int view, int pln, int axis){
	if(mod<15){
		return 24;
	}
	else if(mod==15){
		if(axis==0){
			if(view==0){
				if(pln%3==0)return loli_chnum;
				else        return loli_gridchnum;
			}
			else if(view==1){
				if(pln%3==1)return loli_chnum;
				else        return loli_gridchnum;
			}
            else return 0;
		}
		else if(axis==1){
			if(view==0){
				if(pln%3==1)return loli_plnnum*2;
				else        return loli_plnnum;
			}
			else if(view==1){
				if(pln%3==1)return loli_plnnum*2;
				else        return loli_plnnum;
			}
            else return 0;
		}
        else return 0;
	}
	else if(mod==16){
		return 32;
	}
	else return 0;
}

int INGRID_Dimension::get_plnmax_loli(int mod, int view, int pln, int axis){
	if(mod<15){
		return 11;
	}
	else if(mod==15){
		if(axis==0){
			return loli_plnnum*3;
		}
		else if(axis==1){
			return loli_chnum + loli_gridchnum;
		}
        else return 0;
	}
	else if(mod==16){
		return 18;
	}
    else return 0;
}

bool INGRID_Dimension::get_posi_lolirecon(int mod, int view, int reconpln, int reconch, int axis, double* posxy, double* posz){
	int gridch,grid,pln;
	this->get_plnch_fromrecon_loli(mod, view, reconpln, reconch, axis, &pln, &gridch, &grid);
	this->get_pos_loli_xy(mod, view, pln, gridch, grid, posxy, posz);
	return true;
}

float INGRID_Dimension::get_sciwidth(int mod, int view, int reconpln, int reconch, int axis){
	int gridch,grid,pln;
	this->get_plnch_fromrecon_loli(mod, view, reconpln, reconch, axis, &pln, &gridch, &grid);
	if(axis==0){
		if(grid==0)return loli_scinti_width/2.;  //cm width of scinti
		else return loli_scinti_thick/2.;  //cm thickness of scinti
	}
	else if(axis==1){
		if(grid==0)return loli_scinti_thick/2.;  //cm width of scinti
		else return loli_scinti_width/2.;  //cm thickness of scinti
	}
    else return 0.0;
}

float INGRID_Dimension::get_scithick(int mod, int view, int reconpln, int reconch, int axis){
	int gridch,grid,pln;
	this->get_plnch_fromrecon_loli(mod, view, reconpln, reconch, axis, &pln, &gridch, &grid);
	if(axis==0){
		if(grid==0)return loli_scinti_thick/2.;  //cm width of scinti
		else return loli_scinti_width/2.;  //cm thickness of scinti
	}
	else if(axis==1){
		if(grid==0)return loli_scinti_width/2.;  //cm width of scinti
		else return loli_scinti_thick/2.;  //cm thickness of scinti
	}
    else return 0.0;
}



void INGRID_Dimension::get_loli_gridcell_id(int mod, int view, int pln, int ch, double posx, double posy, double posz,
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


void INGRID_Dimension::get_loli_scintiid_from_cellid(int mod, int view, int pln, int gridcell_id, int* cross_n, int* ch){
	if(view==0){
		if(gridcell_id==0){
			*cross_n=2;
			ch[0] = 0;
			ch[1] = 40;
		}
		else if(gridcell_id>=1 && gridcell_id<=19){
			*cross_n=4;
			ch[0] = 2*gridcell_id -1;
			ch[1] = 2*gridcell_id;
			ch[2] = 40+gridcell_id-1;
			ch[3] = 40+gridcell_id;
		}
		else if(gridcell_id==20){
			*cross_n=2;
			ch[0] = 39;
			ch[1] = 59;
		}
		else if(gridcell_id==21){
			*cross_n=1;
			ch[0] = 60;
		}
		else if(gridcell_id>=22 && gridcell_id<=40){
			*cross_n=2;
			ch[0] = 40+gridcell_id-2;
			ch[1] = 40+gridcell_id-1;
		}
		else if(gridcell_id==41){
			*cross_n=1;
			ch[0] = 79;
		}
	}
	else if(view==1){
		if(gridcell_id==0){
			*cross_n=1;
			ch[0] = 40;
		}
		else if(gridcell_id>=1 && gridcell_id<=19){
			*cross_n=2;
			ch[0] = 40+gridcell_id-1;
			ch[1] = 40+gridcell_id;
		}
		else if(gridcell_id==20){
			*cross_n=1;
			ch[0] = 59;
		}
		else if(gridcell_id==21){
			*cross_n=2;
			ch[0] = 0;
			ch[1] = 60;
		}
		else if(gridcell_id>=22 && gridcell_id<=40){
			*cross_n=4;
			ch[0] = 2*(gridcell_id-21) -1;
			ch[1] = 2*(gridcell_id-21);
			ch[2] = 60+(gridcell_id-21)-1;
			ch[3] = 60+(gridcell_id-21);
		}
		else if(gridcell_id==41){
			*cross_n=2;
			ch[0] = 39;
			ch[1] = 79;
		}
	}
}
#endif
