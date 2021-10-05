#define VME_SLOTS_MAX   21
#define ROCS_MAX        200
#define FADC_CH_NUM     16

#define FADC_TRIG_MAX        1000
#define FADC_HITS_MAX        3000
#define FADC_WAVE_SAMPLE_MAX 100
#define FADC_WAVE_MAX        512

typedef struct
{
  Int_t    fadc_n;
  Int_t    fadc_roc[FADC_TRIG_MAX];
  Int_t    fadc_slot[FADC_TRIG_MAX];
  Double_t fadc_trig_t[FADC_TRIG_MAX];
  Int_t    fadc_trig_n[FADC_TRIG_MAX];
  
  // FADC Waveforms
  Int_t    wave_n;
  Int_t    wave_roc[FADC_WAVE_MAX];
  Int_t    wave_slot[FADC_WAVE_MAX];
  Int_t    wave_ch[FADC_WAVE_MAX];
  Int_t    wave_thr[FADC_WAVE_MAX];
  Float_t  wave_ped[FADC_WAVE_MAX];
  Int_t    wave_nsamples[FADC_WAVE_MAX];
  Int_t    wave[FADC_WAVE_MAX][FADC_WAVE_SAMPLE_MAX];
  
} fadc_evt;

#define TI_N_MAX			500

typedef struct
{
	Int_t		n;
	Int_t		roc[TI_N_MAX];
	Double_t trig_t[TI_N_MAX];
	Int_t		trig_n[TI_N_MAX];
	
	// TI Trigger Type;
  struct
  {
    Int_t n;
    Int_t id[TI_N_MAX];
  } trg;
} ti_evt;

ti_evt tiEvent;
fadc_evt fadcEvent;

typedef struct
{
  const char *name;
  void *address;
  const char *leaf;
} branch_def;

static branch_def evt_branch[] = {

		{"ti_n",      &tiEvent.n,         "ti_n/I"},
		{"ti_roc",    &tiEvent.roc[0],    "ti_roc[ti_n]/I"},
		{"ti_trig_t", &tiEvent.trig_t[0], "ti_trig_t[ti_n]/D"},
		{"ti_trig_n", &tiEvent.trig_n[0], "ti_trig_n[ti_n]/I"},

    {"ti.trg.n",  &tiEvent.trg.n,     "ti.trg.n/I"},
    {"ti.trg.id", &tiEvent.trg.id[0], "ti.trg.id[ti.trg.n]/I"},
    {"fadc_n",        &fadcEvent.fadc_n,            "fadc_n/I"},
    {"fadc_roc",      &fadcEvent.fadc_roc[0],       "fadc_roc[fadc_n]/I"},
    {"fadc_slot",     &fadcEvent.fadc_slot[0],      "fadc_slot[fadc_n]/I"},
    {"fadc_trig_t",   &fadcEvent.fadc_trig_t[0],    "fadc_trig_t[fadc_n]/D"},
    {"fadc_trig_n",   &fadcEvent.fadc_trig_n[0],    "fadc_trig_n[fadc_n]/I"},
       
    {"wave_n",        &fadcEvent.wave_n,            "wave_n/I"},
    {"wave_roc",      &fadcEvent.wave_roc[0],       "wave_roc[wave_n]/I"},
    {"wave_slot",     &fadcEvent.wave_slot[0],      "wave_slot[wave_n]/I"},
    {"wave_ch",       &fadcEvent.wave_ch[0],        "wave_ch[wave_n]/I"},
    {"wave_thr",      &fadcEvent.wave_thr[0],       "wave_thr[wave_n]/I"},
    {"wave_ped",      &fadcEvent.wave_ped[0],       "wave_ped[wave_n]/F"},
    {"wave_nsamples", &fadcEvent.wave_nsamples[0],  "wave_nsamples[wave_n]/I"},
    {"wave",          &fadcEvent.wave[0][0],        Form("wave[wave_n][%d]/I", FADC_WAVE_SAMPLE_MAX)}
	};

void BranchInit(TTree *pT)
{
	// Event data
	for(int i = 0; i < sizeof(evt_branch)/sizeof(evt_branch[0]); i++)
	{
		pT->Branch(evt_branch[i].name, evt_branch[i].address, evt_branch[i].leaf);
		pT->SetBranchAddress(evt_branch[i].name, evt_branch[i].address);
	}
	
//  memset(&vtpCfg, 0, sizeof(vtpCfg));
//  memset(&fadcCfg, 0, sizeof(fadcCfg));
}

//void draw_wave(int entry, int roc, int slot, int ch)
void anafadc_marocsum(int vme_run,int nchip=1,int thisevent=-1)
{

	gStyle->SetOptStat(1);
//		gStyle->SetOptFit(111111);	
	
  TFile *f = new TFile( Form("data/fevme1_%06d.root",vme_run), "READ" );
  TTree *pT = (TTree *)f->Get("T");
//  TTree *pTCfg = (TTree *)f->Get("TCfg");
  //BranchInit(pT, pTCfg);
  BranchInit(pT);

  int tw=48;
//  int tw=96;
  TH1I *pH[14];
  pH[0] = new TH1I("h1","fadc0_empty;ns;adc*Nevent",tw, 0,4*tw);
  pH[1] = new TH1I("h2","fadc1_laser1;ns;adc*Nevent",tw, 0,4*tw);
  pH[2] = new TH1I("h3","fadc2_laser2;ns;adc*Nevent",tw, 0,4*tw);
// pH[3] = new TH1I("h4","fadc3_sum64;ns;adc*Nevent",tw, 0,4*tw);
// pH[4] = new TH1I("h5","fadc4_dynode;ns;adc*Nevent",tw, 0,4*tw);
  pH[3] = new TH1I("h4","fadc3_empty;ns;adc*Nevent",tw, 0,4*tw);
  pH[4] = new TH1I("h5","fadc4_chip2_sum64;ns;adc*Nevent",tw, 0,4*tw);
  pH[5] = new TH1I("h6","fadc5_chip2_quad1;ns;adc*Nevent",tw, 0,4*tw);
  pH[6] = new TH1I("h7","fadc6_chip2_quad2;ns;adc*Nevent",tw, 0,4*tw);
  pH[7] = new TH1I("h8","fadc7_chip2_quad3;ns;adc*Nevent",tw, 0,4*tw);
  pH[8] = new TH1I("h9","fadc8_chip2_quad4;ns;adc*Nevent",tw, 0,4*tw);
  pH[9] = new TH1I("h10","fadc9_chip1_sum64;ns;adc*Nevent",tw, 0,4*tw);
  pH[10] = new TH1I("h11","fadc10_chip1_quad1;ns;adc*Nevent",tw, 0,4*tw);
  pH[11] = new TH1I("h12","fadc11_chip1_quad2;ns;adc*Nevent",tw, 0,4*tw);
  pH[12] = new TH1I("h13","fadc12_chip1_quad3;ns;adc*Nevent",tw, 0,4*tw);
  pH[13] = new TH1I("h14","fadc13_chip1_quad4;ns;adc*Nevent",tw, 0,4*tw);
  
  Long64_t nevent = pT->GetEntries();
  printf("pT->GetEntries() = %lld\n", nevent);

  int bad=0;
  for(Long64_t entry=2;entry<nevent-11;entry++) {
     pT->GetEntry(entry);
	for(int i=5;i<9;i++){
         for(int j=0;j<fadcEvent.wave_nsamples[i];j++) {
//		cout << 4095-fadcEvent.wave[i][j] << endl;
//                if(4095-fadcEvent.wave[i][j]<430) {bad++; cout<< "bad " << entry << " " <<  i << " " << j << endl;}
         }	
	}
  }
  cout << bad << " bad events" << endl;

//  for(Long64_t entry=2;entry<nevent-11;entry++)  // 1st 2 and last 11 events are junk
//  for(Long64_t entry=2;entry<nevent-2;entry++)  // 1st 2 and last 2 events are junk
//  for(Long64_t entry=2;entry<1e5+2;entry++)  // 1st 2 and other events are junk
//  for(Long64_t entry=10;entry<11;entry++)
  int event1,event2;
  if (thisevent==-1) {event1=2;event2=nevent-2;}
  else  {event1=thisevent;event2=thisevent+1;}
  for(Long64_t entry=event1;entry<event2;entry++)
  {
//    if(!(entry%1000)) printf("%lld ", entry);

  pT->GetEntry(entry);
//  cout << fadcEvent.wave_n << endl;

      for (int i=0;i<fadcEvent.wave_n;i++){
//    if(fadcEvent.wave_roc[i] != roc) continue;
//    if(fadcEvent.wave_slot[i] != slot) continue;
//    if(fadcEvent.wave_ch[i] != ch) continue;

//    for(int j=0;j<fadcEvent.wave_nsamples[i];j++)  pH[i]->SetBinContent(j+1,fadcEvent.wave[i][j]);
//	cout << fadcEvent.wave_n << endl;
//	cout << fadcEvent.wave_nsamples[i] << endl;

/*
cable connection after 2020/07/07, when testing maroc sum production board
fadc 0, empty
fadc 1, laser signal at 2V
fadc 2, laser signal with delay at 0.5V
fadc 3, empty    (green cable)
fadc 4, chip 2, maroc sum64  (black cable) 
fadc 5, chip 2, maroc Quad1  (white cable 03)
fadc 6, chip 2, maroc Quad2  (white cable 04)
fadc 7, chip 2, maroc Quad3  (white cable 05)
fadc 8, chip 2, maroc Quad4  (white cable 06)
fadc 9, chip 1, maroc sum64 
fadc 10, chip 1, maroc Quad1
fadc 11, chip 1, maroc Quad2
fadc 12, chip 1, maroc Quad3
fadc 13, chip 1, maroc Quad4
and 
maroc sum64 signal is postive, all quad are negative
*/

    for(int j=0;j<fadcEvent.wave_nsamples[i];j++) {
//	 if (i<3) pH[i]->Fill(4*j,4095-fadcEvent.wave[i][j]);
//	 else pH[i]->Fill(4*j,fadcEvent.wave[i][j]);

	if (i==0 || i==1 || i==2 || i==3 || i==5 || i== 6 || i==7 || i== 8|| i==10 || i== 11|| i==12 || i== 13) pH[i]->Fill(4*j,4095-fadcEvent.wave[i][j]);
	else pH[i]->Fill(4*j,fadcEvent.wave[i][j]);
	pH[i]->SetBinError(j+1,0);
    }
  }

  }


  TCanvas *tC = new TCanvas("c2","c2",1200,1000);
  tC->Divide(1,2,0.01,0.01);

    tC->cd(1);
    pH[1]->Draw();
    tC->cd(2);
    pH[2]->Draw();


  TCanvas *pC = new TCanvas("c1","c1",1200,500);
//  pC->Divide(1,3,0.01,0.01);

    pC->cd(1);
    pH[0]->SetMinimum(0);
//    pH[0]->SetMaximum(4095);
    pH[0]->Draw("Bar0");
//pC->SaveAs(Form("data/fevme1_%06d.png",vme_run));

  TCanvas *simpleC = new TCanvas("c3","c3",1600,1000);
  if (nchip==1){
	 simpleC->Divide(2,3,0.01,0.01);
	  for(int i=4;i<9;i++)
	  {
	    simpleC->cd(i-2);
	// 	if (i>4) pH[i]->SetMinimum(2300);
	//	else     pH[i]->SetMinimum(1400);
	//	if (i>4) pH[i]->SetMaximum(3000*1e5);
	//	else     pH[i]->SetMaximum(3000*1e5);

	    pH[i]->SetMinimum(0);
	//    pH[0]->SetMaximum(4095);
	    pH[i]->Draw("Bar0");
	  }
  }
  else if (nchip==2) {
	simpleC->Divide(5,2,0.01,0.01);

         for(int i=4;i<14;i++)
          {
            simpleC->cd(i-3);

            pH[i]->SetMinimum(0);
        //    pH[0]->SetMaximum(4095);
            pH[i]->Draw("Bar0");

	 }	
  }
  else cout << "not sue how many chips";

simpleC->SaveAs(Form("data/fevme1_%06d.png",vme_run));

}

