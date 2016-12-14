#include"Midas_to_Root_new.hxx"


FileStat_t                                  fs;
INGRID_Ch_config*            fINGRID_Ch_config;
INGRID_BadCh_mapping*    fINGRID_BadCh_mapping;
INGRID_Dimension*            fINGRID_Dimension;
//_____________________________________________

//____________________________________________
Bool_t flana;
void Analysis(vector<Hit> hit[][cNumCyc]);
void Book();
void Read(ND::TND280RawEvent* re, vector<Hit> hit[][cNumCyc]);
void ProcessFile(const char *FileName);

//_______________________________________________________________
int main(int argc,char *argv[]){
  TROOT root("GUI","GUI");
  TApplication theApp("App",0,0);
  int c=-1;
  char FileName[300],root_file_name[300];
  int run_number, calib_run_number;
  cAnaTrg    =  cTrgBeam; //default Beam trigger
  cAnaEvt    =       100;
  cAnaMod    =        -1;
  cAnaCyc    =        -1;
  while ((c = getopt(argc, argv, "hr:c:t:b:m:l:")) != -1) {
    switch(c){
    case 'r':
      run_number=atoi(optarg);
      //sprintf(FileName,"%s/ingrid_%08d_0000.daq.mid",data_file_folder,run_number);
      sprintf(FileName,"%s/ingrid_%08d_0000.daq.mid.gz",data_file_folder,run_number);
      break;
    case 'c':
      calib_run_number=atoi(optarg);
      break;
    case 't':
      cAnaTrg=atoi(optarg);
      break;
    case 'b':
      cAnaEvt=atoi(optarg);
      break;
    case 'l':
      cAnaCyc=atoi(optarg);
      break;
    case 'h':
      cout<<"-r [run number]\t midas file to root file"<<endl;
      break;
    case 'm':
      cAnaMod=atoi(optarg);
      break;
    case '?':
      cout<<"Unknown option"<<endl;
      cout<<"-r [run number]\t midas file to root file"<<endl;
      exit(1);
      break;
    }
  }//option end
  fINGRID_Ch_config     = new INGRID_Ch_config();
  fINGRID_BadCh_mapping = new INGRID_BadCh_mapping();
  fINGRID_Dimension     = new INGRID_Dimension();
  read_calib(calib_run_number);
  Book(run_number);
  ProcessFile(FileName);
  Write();
  return 0;
}

//_________________________________________________________________

void ProcessFile(const char *FileName) {
  if ( gSystem->GetPathInfo(FileName,fs) ) {
    std::cerr << "Cannot find file: " << FileName << std::endl;
    return;
  }
  ND::TMidasFile mf;
  mf.Open(FileName);
  cout<<"loop all event..."<<endl;

  while ( ND::TND280RawEvent* re = mf.ReadRawEvent()) {
    re->PromoteMidasBanks(false);

    for(Int_t i=0;i<cNumMod;i++){
      for(Int_t j=0;j<cNumCyc;j++){
	hit[i][j].clear();
      }
    }
    TrgId=-1;

    if(re->size()!=0){//eliminate header and laster
      Read(re, hit);
      if(TrgId==cAnaTrg)Analysis(hit);
 
    }

    delete re;
    if(NumEvt==cAnaEvt)break;
    if(NumEvt%100==0)cout<<"event:"<<NumEvt<<endl;
    if(NumEvt%10000==0)Save();
  }// End loop over events
}
//_________________________________________________________________
#define DEBUG 0
void Read(ND::TND280RawEvent* re, vector<Hit> hit[][cNumCyc]) {
  //Header
  ND::TRawDataHeader header =re->GetHeader();
  UTime   = header.GetTimeStamp();

  //TRunInfo
  ND::THandle<ND::TRunInfoBank> RunInfoBank;
  while ( RunInfoBank = re->GetMidasBank<ND::TRunInfoBank>("XRUN",RunInfoBank) ) {
    ND::TRunInfoBank& runinfo = re->UseMidasBank<ND::TRunInfoBank>("XRUN");
    NumEvt = runinfo.GetSeqNumber();
  }

  //MCMBank
  ND::THandle<ND::TMCMBank> mcmBank;
  while ( mcmBank = re->GetMidasBank<ND::TMCMBank>("IMCM",mcmBank) ) {
    ND::TMCMBank& mcm = re->UseMidasBank<ND::TMCMBank>("IMCM");
    fTrgTime = mcm.GetUnixTimeSSecTrig();

    //### for study MCM time problem ######  
    MCMTime_a_now        =   mcm.GetUnixTimeTrig()  * constant;
    MCMTime_b_now        =   mcm.GetUnixTimeSSecTrig();
    if(NumEvt <= 1  ){
      MCMTime_a_before   =   MCMTime_a_now;  
      MCMTime_b_before   =   MCMTime_b_now;  
    }
    if(NumEvt >  1 ){
      MCMTime_c          =   mcm.GetTriggerSep();
      //Print();
      MCMTime_diff       =   ( MCMTime_a_now + MCMTime_b_now ) - (MCMTime_a_before + MCMTime_b_before );
      MCMTime_a_before   =   MCMTime_a_now;  
      MCMTime_b_before   =   MCMTime_b_now;  
    }
 
  }
  //Trigger Bank
  ND::THandle<ND::TTriggerBank> triggerBank;
  while ( triggerBank = re->GetMidasBank<ND::TTriggerBank>("ITRI",triggerBank) ) {
    ND::TTriggerBank& trigger = re->UseMidasBank<ND::TTriggerBank>("ITRI");
    TrgId  = (trigger.GetTriggerWord()>>48) & 0xffff;
    nSpill = (trigger.GetTriggerWord()>>32) & 0xffff;
  }

  //if(DEBUG&&nSpill>16915)cout<<"------TrgTime:"<<fTrgTime*10<<" --------"<<endl; 
 
  if(TrgId==cAnaTrg){
    ND::THandle<ND::TTripTDigitBank> triptBank;
    while ( triptBank = re->GetMidasBank<ND::TTripTDigitBank>("",triptBank) ) {
      // Create an iterator over digits
      
      ND::TMidasTripTDigitItr itr(triptBank->GetMidasTripTDigitItr());
      while ( ! itr.EOD() ) {
	ND::TMidasTripTDigit digit(itr.Get());
	Int_t rmm     =  digit.GetRMMNum();
	Int_t tfb     =  digit.GetTFBNum();
	Int_t trip    =  digit.GetTripTNum();
	Int_t trip_ch =  digit.GetChannelNum();
	Int_t cycle   =  digit.GetIntegrationNum();
	Int_t mod,view,plane,ch;
	if(fINGRID_Ch_config->channel_configuration(&rmm,&tfb,&trip,&trip_ch,&mod,&view,&plane,&ch)&&(cAnaMod==-1||cAnaMod==mod)){
	  long tdc = digit.GetTimeOffset();
	  bool tpl,veto;
	  int tmp;int tmpch;
	  if(plane>10){tpl=false;veto=true;tmp=plane-11;}
	  else {tpl=true;veto=false;tmp=plane;}
	  if(view==0&&tpl){tmpch=ch+24;}
	  else if(view==1&&tpl){tmpch=ch;}
	  else if(veto){tmpch=ch;}
	  if(tdc<16777201&&!fINGRID_BadCh_mapping->badchannel(&mod, &tmp, &tmpch, &tpl, &veto)){
	    const ND::TRawDataHeader& h = re->GetHeader();int evno = h.GetSeqNo();    // Event Sequence Number
	    int iint = digit.GetIntegrationNum();    // = Capacitor number
	    int icoff = triptBank->GetTFBStartCycle();int co = iint - icoff;
	    if (co<0) co += 23;
	    Hit fhit;
	    fhit.pln    = plane;
	    fhit.ch     =    ch;
	    fhit.view   =  view;
	    fhit.adc    = digit.GetHighGainADC();  
	    fhit.tdc    = digit.GetTimeOffset() + digit.GetTimeOffsetT0();
	    fhit.t0     = digit.GetTimeOffsetT0();
	    fhit.rawtdc = digit.GetTimeOffset();
	    //cout<<"RMM:"<<rmm<<" TFB:"<<tfb<<"\tcycle:"<<cycle<<"\tco"<<co<<"\tT0; "<<digit.GetTimeOffsetT0()<<"\ttdc"<<fhit.tdc<<endl;
	    hit[mod][co].push_back(fhit);
	    /*
	    if(TrgId==1||TrgId==2)
	      hit[mod][co].push_back(fhit);
	    else if(TrgId==128)
	      hit[mod][iint].push_back(fhit);
	    */
	  }
	}//Ch_config
      }//itr
    }//tribBank
  }//TrgId

}

//_________________________________________________________________


#define Debug 0
void Analysis(vector<Hit> hit[][cNumCyc]){
  for(Int_t numcyc=0;numcyc<cNumCyc;numcyc++){
    //for(Int_t numcyc=19;numcyc<20;numcyc++){
    for(Int_t nummod=0;nummod<cNumMod;nummod++){if(cAnaMod==-1||cAnaMod==nummod){


      //for(Int_t nummod=6;nummod<7;nummod++){
      if(Debug)cout<<"----nummod: "<<nummod<<"-------"<<endl;
      fSortTdc(hit[nummod][numcyc]);
      tdc_calib(hit, nummod, numcyc);
      adc_calib(hit, nummod, numcyc);
      

      vector<Hit> hitclster;
      Int_t nClster=0;
      Long_t time;

      adc.  clear();
      pe.   clear();
      nsec. clear();
      tdc.  clear();
      view. clear();
      pln.  clear();
      ch.   clear();
      nHit = hit[nummod][numcyc].size();
      for(Int_t i=0;i<nHit;i++){
	adc. push_back(hit[nummod][numcyc][i].adc );
	pe.  push_back(hit[nummod][numcyc][i].pe  );
	tdc. push_back(hit[nummod][numcyc][i].tdc );
	nsec.push_back(hit[nummod][numcyc][i].time);
	pln. push_back(hit[nummod][numcyc][i].pln );
	view.push_back(hit[nummod][numcyc][i].view);
	ch.  push_back(hit[nummod][numcyc][i].ch  );
	
      }
      Cycle    = numcyc;
      fMod     = nummod;
      tree->Fill();
      

      //if(NumEvt==57){cout<<"Mod:"<<fMod<<" hit"<<nHit<<endl;Print(hit,nummod,numcyc);}

      /*
      while(fFindTimeClster(hit,hitclster,nummod,numcyc,time)){//hit -> hit cluster
	fSortPln(hitclster);
	adc_calibcls(hitclster, nummod, numcyc);
	//Fill
	pe.   clear();
	nsec. clear();
	view. clear();
	pln.  clear();
	ch.   clear();
	nHit = hitclster.size();
	for(Int_t i=0;i<nHit;i++){
	  pe.  push_back(hitclster[i].pe  );
	  nsec.push_back(hitclster[i].time);
	  pln. push_back(hitclster[i].pln );
	  view.push_back(hitclster[i].view);
	  ch.  push_back(hitclster[i].ch  );
	}
	Cycle    = numcyc;
	fTime    = time;
	fMod     = nummod;
	tree->Fill();
	nClster++;
	hitclster.clear();
      }//clustering
      */
      }}//nummod
    }//numcyc
}//Analysis


