#include <iostream>
#include <TString.h>
#include <fstream>
#include <TChain.h>
#include <TMVA/Factory.h>
#include <TMVA/Tools.h>
using std::cout;
using std::endl;

TCut mycuts, mycutb;

void config( TMVA::Factory * factory, TString id )
{
   cout << "\e[31mbegin to read config\e[0m" <<endl;

   //TCut mycuts, mycutb; 
   int nSignalTrain = 0, nBckTrain = 0;

   TChain * signalTrainChain = new TChain("event");
   TChain * bckgndTrainChain = new TChain("event");

   TChain * signalTestChain = new TChain("event");
   TChain * bckgndTestChain = new TChain("event");

   TChain * signalWeightTrainChain = new TChain("event");
   TChain * signalWeightTestChain = new TChain("event");

   TChain * bckgndChain = new TChain("event");

   std::ifstream cfg( id+".conf" );

   while( cfg.good() ){
      TString key, vals;
      cfg >> key;
      vals.ReadLine( cfg, kFALSE );
      std::vector<TString> tokens = TMVA::Tools::Instance().SplitString(vals, ',');

      cout << key << ": " << vals <<endl;

      if( key[0] == '#' )continue;

      if( key == "variable"){
         cout << "variables: " << endl;
         for( size_t i=0; i < tokens.size(); ++i ){
            TString var = tokens[i].Strip(TString::kBoth, ' ');
            if( var.IsNull() )continue;
            factory->AddVariable( var, 'F' );
            cout << var << endl;
         }
      }
      if( key == "spectator" ){
         cout << "spectators: " << endl;
         for( size_t i=0; i < tokens.size(); ++i ){
            TString var = tokens[i].Strip(TString::kBoth, ' ');
            if( var.IsNull() )continue;
            factory->AddSpectator( var, 'F' );
            cout << var << endl;
         }
      }
      if( key == "cut" )
         mycuts = mycutb = vals;
      if( key == "cut_signal" )
         mycuts = vals;
      if( key == "cut_background" )
         mycutb = vals;
      if( key == "weight" )
         factory->SetWeightExpression( vals, "" );

      if( key == "weightS" )
         factory->SetSignalWeightExpression(vals); 
      if( key == "weightB" )
         factory->SetBackgroundWeightExpression(vals);

      if( key == "signal_trainfile")
         signalTrainChain->Add( vals.Remove(TString::kBoth,' ') );
      if( key == "signalWight_trainfile"){
         signalWeightTrainChain->Add( vals.Remove(TString::kBoth,' ') );
         signalTrainChain->AddFriend(signalWeightTrainChain);    
      }
      if( key == "background_trainfile" )
         bckgndTrainChain->Add( vals.Remove(TString::kBoth,' ') );

      if( key == "signal_train_n" )
         nSignalTrain = vals.Atoi();
      if( key == "background_train_n" )
         nBckTrain = vals.Atoi();
      if( key == "background_file" )
         bckgndChain->Add( vals.Remove(TString::kBoth,' ') );
 

      
      if( key == "signal_testfile")
         signalTestChain->Add( vals.Remove(TString::kBoth,' ') );
      if( key == "background_testfile" )
         bckgndTestChain->Add( vals.Remove(TString::kBoth,' ') );

      if( key == "signalWight_testfile"){
         signalWeightTestChain->Add( vals.Remove(TString::kBoth,' ') );
         signalTestChain->AddFriend(signalWeightTestChain);    
      }

   }
   

   cout << "signal taining events: " << signalTrainChain->GetEntries() <<endl;
   //cout << "backgr training events: " << bckgndTrainChain->GetEntries() <<endl;
   cout << "signal test events: " << signalTestChain->GetEntries() <<endl;
   //cout << "backgr test events: " << bckgndTestChain->GetEntries() <<endl;
   cout << "backgr events: " << bckgndChain->GetEntries() <<endl;

   //if( signalTrainChain->GetEntries() * bckgndTrainChain->GetEntries() == 0 )return;
   //if( signalTestChain->GetEntries() * bckgndTestChain->GetEntries() == 0 )return;
  if( signalTrainChain->GetEntries() * bckgndChain->GetEntries() * signalTestChain->GetEntries()  == 0 )return;

cout<<"=========================TaoJQ================"<<endl;
  factory->AddSignalTree( signalTrainChain, 1.0,  "train");
  //factory->AddBackgroundTree( bckgndTrainChain, 1.0, "train" );
  //factory->AddBackgroundTree(bckgndChain, 1.0, "train" );

  factory->AddSignalTree( signalTestChain, 1.0,  "test");
  //factory->AddBackgroundTree( bckgndTestChain, 1.0, "test" );
  //factory->AddBackgroundTree(bckgndChain, 1.0, "test" );

  factory->AddBackgroundTree( bckgndChain, 1.0 );


/*
   if( mycuts == mycutb )
      factory->PrepareTrainingAndTestTree(
         mycuts,
         nSignalTrain,nBckTrain,
         0,0,
         "SplitMode=Random:NormMode=NumEvents:!V"
         );
   else
      factory->PrepareTrainingAndTestTree(
            mycuts,
            mycutb,
            //"nTrain_Signal=100000:nTrain_Background=100000:nTest_Signal=100000:nTest_Background=100000:SplitMode=Random:NormMode=NumEvents:!V"  
            "SplitMode=Random:NormMode=NumEvents:!V" /// TODO add Num of training support /
            );
*/
}

vector<std::string> getVaribles( TString id )
{
   vector<std::string> res;
   std::ifstream cfg( id+".conf" );

   while( cfg.good() ){
      TString key, vals;
      cfg >> key;
      vals.ReadLine( cfg, kFALSE );

      //clog << key << endl;
      if( key=="variable" ){
         std::vector<TString> tokens = TMVA::Tools::Instance().SplitString(vals, ',');
         for( size_t i=0; i < tokens.size(); ++i ){
            TString var = tokens[i].Strip(TString::kBoth, ' ');
            if( var.IsNull() )continue;
            res.push_back( std::string(var.Data()) );
            clog << "var"<<i <<": " << var << ' ';
         }
         clog << '\n';
         return res;
      }
   }
   return res;
}

