//
//###########################################################
//###            MIDAS file -> INGRID ROOT File           ###
//###########################################################
//###                Made by M.Otani                      ###
//###########################################################
//#include"DSTMaker.hxx"

#include"DSTMaker_Optimized.hxx"
void Read       (ND::TND280RawEvent* re);
void ProcessFile(const char *FileName);

int run_number, subrun;
//_______________________________________________________________
int main(int argc,char *argv[]){
  TROOT root("GUI","GUI");
  TApplication theApp("App",0,0);
  int c=-1;
  char FileName[300]      ="hoge.mid";  //Midas file name
  //char root_file_name[300]="hoge.root"; //output ROOT file name
  fINGRID_Ch_config       = new INGRID_Ch_config();
  fPM_Ch_config           = new PM_Ch_config();
  fINGRID_BadCh_mapping   = new INGRID_BadCh_mapping();
  fINGRID_Dimension       = new INGRID_Dimension();
  cosmic     =  false;  
  run_number =  0; 
  subrun     =  0;          //ingrid_${run_number}_{subrun}.root
  cAnaEvt    = -1;          //default, all event will be processed
  PModule    =  false;
  while ((c = getopt(argc, argv, "hr:s:t:b:p")) != -1) {
    switch(c){
    case 'r':
      run_number=atoi(optarg);
      break;
    case 's':
      subrun=atoi(optarg);
      break;
    case 't':
      cAnaTrg=atoi(optarg);
      break;
    case 'b':
      cAnaEvt=atoi(optarg);
      break;
    case 'p':
      PModule = true;
      break;
    case 'h':
      cout<<"-r [run number]\t midas file to root file"<<endl;
      break;
    case '?':
      cout<<"Unknown option"<<endl;
      cout<<"-r [run number]\t midas file to root file"<<endl;
      exit(1);
      break;
    }
  }//option end
  if(cAnaTrg==128)
    cosmic = true;

  //sprintf(FileName,"/export/scraid2/data/ingrid/daqdata/ingrid_%08d_%04d.daq.mid.gz",run_number,subrun);
  //sprintf(FileName,"/export/scraid2/data/bquilain/Midas_All/ingrid_%08d_%04d.daq.mid.gz",run_number,subrun);
  sprintf(FileName,"/home/bquilain/ingrid_%08d_%04d.daq.mid.gz",run_number,subrun);
  //sprintf(FileName,"%s/ingrid_%08d_%04d.daq.mid.gz",data_file_folder,run_number,subrun);
  if ( gSystem->GetPathInfo(FileName,fs) ) {
    std::cerr << "Cannot find file: " << FileName << std::endl;
    /*
    sprintf( FileName, "/data/daqdata/ingrid_%08d_%04d.daq.mid.gz",run_number,subrun);
    if ( gSystem->GetPathInfo(FileName,fs) ) {
      std::cerr << "Cannot find file: " << FileName << std::endl;
      
    }
    */
    exit(1);
  }
  Book(run_number, subrun);

  ProcessFile(FileName);
  Write();
  return 0;
}


//_________________________________________________________________


