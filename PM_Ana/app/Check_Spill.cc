#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<sys/stat.h>
using namespace std;

int main(int argc, char *argv[]){
  struct stat st;
  int             c  =     -1;

  char IngFile[300], BSDFile[300];
  while ((c = getopt(argc, argv, "i:b:")) != -1) {
    switch(c){
    case 'i':
      sprintf( IngFile, "%s", optarg);
      break;
    case 'b':
      sprintf( BSDFile, "%s", optarg);
      break;
    }
  }

  if( stat(IngFile, &st) ){
    cout << "not exist " << IngFile << endl;
    exit(1);
  }

  if(  stat(BSDFile, &st) ){
    cout << "not exist " << BSDFile << endl;
    exit(1);
  }


  vector<pair<int,int> > ing_spill_list;
  vector<pair<int,int> > bsd_spill_list;
  vector<pair<long,long> > ing_nrun_list;
  vector<pair<int,int> > bsd_nrun_list;
  ing_spill_list.clear();
  bsd_spill_list.clear();
  ing_nrun_list.clear();
  bsd_nrun_list.clear();

  long t1, t;
  ifstream fBSDFile(BSDFile);
  while(fBSDFile>>t>>t1){
    pair<int,int> t2, t3;
    t2.first  = 0;
    t2.second = t1;
    t3.first  = 0;
    t3.second = t;
    bsd_spill_list.push_back(t2);
    bsd_nrun_list.push_back(t3);
  }
  fBSDFile.close();
  cout << "read... " << bsd_spill_list.size();

  ifstream fIngFile(IngFile);
  while(fIngFile>>t>>t1){
    pair<long,long> t2,t3;
    t2.first  = 0;
    t2.second = t1;
    t3.first  = 0;
    t3.second = t;
    ing_spill_list.push_back(t2);
    ing_nrun_list.push_back(t3);
  }
  fIngFile.close();
  cout << "read... " << ing_spill_list.size();

  int lastj=0;
  for(int i=0; i<bsd_spill_list.size(); i++){
    if(i%10000==0)cout << "check :"<< i << endl;
    for(int j=lastj; j<ing_spill_list.size(); j++){
      if( bsd_spill_list[i].second == ing_spill_list[j].second ){
	bsd_spill_list[i].first = 1;
	ing_spill_list[j].first = 1;
	bsd_nrun_list[i].first = 1;
	ing_nrun_list[j].first = 1;
	lastj = j;
	break;
      }//find
    }//ing loop
  }

  int exist_ing_nexist_bsd=0;
  for(int i=0; i<ing_spill_list.size(); i++){
    if( ing_spill_list[i].first == 0 ){
      exist_ing_nexist_bsd++;
      cout << ing_spill_list[i].second << endl;
    }
  }
  ofstream wfile("MissSpill.txt");
  int nexist_ing_exist_bsd=0;
  for(int i=0; i<bsd_spill_list.size(); i++){
    if( bsd_spill_list[i].first == 0 ){
      nexist_ing_exist_bsd++;

      wfile << bsd_nrun_list[i].second  << "\t"
	    << bsd_spill_list[i].second << endl;
    }
  }
  wfile.close();

  cout << "Good Spill at Ingrid     :" << ing_spill_list.size() << endl;
  cout << "Good Spill at BSD        :" << bsd_spill_list.size() << endl;
  cout << "Exist BSD but Not Ingrid :" << nexist_ing_exist_bsd  << endl;
  cout << "Exist Ingrid but Not BSD :" << exist_ing_nexist_bsd  << endl;
  cout << "create MissSpill.txt" << endl;

}
