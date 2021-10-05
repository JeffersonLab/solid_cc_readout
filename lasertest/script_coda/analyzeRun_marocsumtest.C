
void analyzeRun_marocsumtest(int vme_run) {
	
	gStyle->SetOptStat(111111);
		gStyle->SetOptFit(111111);	
	
//	int vme_run = 145;
	
//before 10x amp added on 2019/10/15	
//	int nped_a = 12;
//	int nped_b = 16;
//	int nped   =  nped_b-nped_a;              // number of samples used for pedestal determination
//	
//	// peak around bin 22, full width about 6 bin
//	int nint_a = 16;
//	int nint_b = 30;
//	int nint   = nint_b - nint_a; // number of samples used for integral determination
	
//after 10x amp added on 2019/10/15, use integral width as multiple of 3 to reduce noise of 3 bin wide
//	int nped_a = 11;
//	int nped_b = 20;
//	int nped   =  nped_b-nped_a;              // number of samples used for pedestal determination
//	
//	// peak around bin 27, full width about 6 bin
//	int nint_a = 20;
//	int nint_b = 35;
//	int nint   = nint_b - nint_a; // number of samples used for integral determination
	

//after 2019/10/24, add 450ohm resistor in series with MAROC sum output.this,along with FADC 50ohm termination will divide signal by 10, ~10x NIM amp is bypassed
//use integral width as multiple of 3 to reduce noise of 3 bin wide
//	int nped_a = 7;
//	int nped_b = 16;
//	int nped   =  nped_b-nped_a;              // number of samples used for pedestal determination
//	
//	// peak around bin 22, full width about 6 bin
//	int nint_a = 16;
//	int nint_b = 31;
//	int nint   = nint_b - nint_a; // number of samples used for integral determination

//after 2019/10/25, add onboard 0.5x amp and impedance, sum of 64 with signal reversed which can be subtracted from 4095 
//use integral width as multiple of 3 to reduce noise of 3 bin wide
	int nped_a = 13;
	int nped_b = 22;
	int nped   =  nped_b-nped_a;              // number of samples used for pedestal determination
	
	// peak around bin 96, full width about 6 bin
	int nint_a = 22;
	int nint_b = 37;
	int nint   = nint_b - nint_a; // number of samples used for integral determination


	TFile *fIn = new TFile( Form("data/fevme1_%06d.root",vme_run), "READ" );
	TTree *T   = (TTree*)fIn->Get( "T" );
	
	int nwave;
	int wave[10][100];
	T->SetBranchAddress( "wave_n", &nwave );
	T->SetBranchAddress( "wave"  , &wave  );
	
	
	
	TH1I *hped = new TH1I( "hped", "Pedestal",   4095, 0., 4095. );
	hped->SetTitle( Form("Pedestal Distribution (average fADC count between %d - %d ns)", 
		nped_a*4, nped_b*4) );
	hped->GetXaxis()->SetTitle( "Pedestal [fADC counts]" );
	
	TH1I *hamp    = new TH1I( "hamp",    "Amplitude",                       4095, 0., 4095. );
	TH1I *hampSub = new TH1I( "hampSub", "Amplitude (pedestal subtracted)", 4500, -500., 4000. );
	hamp->SetTitle( "Amplitude Distribution (no pedestal subtraction)" );
	hamp->GetXaxis()->SetTitle( "Amplitude [fADC Counts]" );
	hampSub->SetTitle( "Amplitude Distribution (with pedestal subtraction)" );
	hampSub->GetXaxis()->SetTitle( "Amplitude [fADC]" );
	
	
	TH1I *hint    = new TH1I( "hint",    "Integral",                       40500, -500., 40000. );
	TH1I *hintSub = new TH1I( "hintSub", "Integral (pedestal subtracted)", 40500, -500., 40000. );
	hint->SetTitle( "Integral Distribution (no pedestal subtraction)" );
	hint->GetXaxis()->SetTitle( "Integral [fADC Counts]" );
	hintSub->SetTitle( "Integral Distribution (with pedestal subtraction)" );
	hintSub->GetXaxis()->SetTitle( "Integral [fADC]" );
	
	
	cout << "N of events " << T->GetEntries() << endl;
	for( int ievt = 0; ievt < T->GetEntries(); ievt++ ) {
	  
	  
	  if( ievt%100000 == 0 ) cout << "Processing event " << ievt << endl;
	  T->GetEvent(ievt);
	  
	  
	  
	  int loc_ped = 0;
	  for( int isamp = nped_a; isamp < nped_b; isamp++ ) loc_ped += 4095-wave[0][isamp];
	  
	  double loc_amp = 0.; double loc_int = 0.;
	  for( int isamp = nint_a; isamp < nint_b; isamp++ ) {	    
	    double f_wave = static_cast<double>( 4095-wave[0][isamp] );
	    loc_int += f_wave;
	    if( f_wave > loc_amp ) loc_amp = f_wave;
	  }
	  
	  double f_ped = static_cast<double>( loc_ped ) / static_cast<double>( nped );
	  
	  hped->Fill( f_ped );
	  hamp->Fill( loc_amp );
	  hint->Fill( loc_int );
	  
	  loc_amp = loc_amp - f_ped;
	  loc_int = loc_int - f_ped*static_cast<double>( nint );
	  
	  hampSub->Fill( loc_amp );
	  hintSub->Fill( loc_int );
	  
	}
	
	
	hamp->Rebin(2);
	hampSub->Rebin(2);
	hint->Rebin(4);
	hintSub->Rebin(4);
	
	hped->SetLineWidth( 2 );
	hamp->SetLineWidth( 2 );
	hampSub->SetLineWidth( 2 );
	hint->SetLineWidth( 2 );
	hintSub->SetLineWidth( 2 );
	
	hped->SetLineColor( kBlack );
	hamp->SetLineColor( kBlack );
	hampSub->SetLineColor( kBlack );
	hint->SetLineColor( kBlack );
	hintSub->SetLineColor( kBlack );
	
	TCanvas *c1 = new TCanvas( "c1", "ped", 1000, 1000 );
	c1->SetTickx(); c1->SetTicky();
	c1->cd(); hped->Draw();
	
	TCanvas *c2 = new TCanvas( "c2", "amp", 1000, 1000 );
	c2->SetTickx(); c2->SetTicky();
	c2->cd(); hamp->Draw();
	
	TCanvas *c3 = new TCanvas( "c3", "int", 1000, 1000 );
	c3->SetTickx(); c3->SetTicky();
	c3->cd(); hint->Draw();
	
	TCanvas *c4 = new TCanvas( "c4", "ampSub", 1000, 1000 );
	c4->SetTickx(); c4->SetTicky();
	c4->cd(); hampSub->Draw();
//	gPad->SetLogy(1);
	
	TCanvas *c5 = new TCanvas( "c5", "intSub", 1000, 1000 );
	c5->SetTickx(); c5->SetTicky();
	c5->cd(); hintSub->Draw();
//	gPad->SetLogy(1);
	
	
	  TF1 *fPed = new TF1( "fPed", "gaus", 0., 1. );

	  double cmax, xmax;	  

	  hampSub->GetXaxis()->SetRangeUser( 20., 4000. );
	  cmax = static_cast<double>( hampSub->GetMaximum() );
	  xmax = static_cast<double>( hampSub->GetBinCenter(hampSub->GetMaximumBin()) );
  	  fPed->SetRange( xmax-(xmax-100), xmax+(xmax-100) );
	  fPed->SetParameters( cmax, xmax, (xmax-100)/2. );
	  hampSub->Fit( "fPed", "R" );

	  hintSub->GetXaxis()->SetRangeUser( 200., 10000. );
	  cmax = static_cast<double>( hintSub->GetMaximum() );
	  xmax = static_cast<double>( hintSub->GetBinCenter(hintSub->GetMaximumBin()) );
	  fPed->SetRange( xmax-(xmax-200), xmax+(xmax-200) );
	  fPed->SetParameters( cmax, xmax, (xmax-200)/2. );	  
	  hintSub->Fit( "fPed", "R" );
	  
//	  TF1 *locF = (TF1*)fPed->Clone( Form("fPed_%d",ich) );
//	  locF->SetLineColor(kRed);
//	  locF->Draw("same");
	
//	hped->GetXaxis()->SetRangeUser(     500.,   700. );
//	hamp->GetXaxis()->SetRangeUser(     500.,   3000. );
//	hampSub->GetXaxis()->SetRangeUser(  -50.,   2500. );
//	hint->GetXaxis()->SetRangeUser(    8000.,  20000. );
//	hintSub->GetXaxis()->SetRangeUser( -500.,  10000. );

	hped->GetXaxis()->SetRangeUser(     0.,   1500. );
	hamp->GetXaxis()->SetRangeUser(     0.,   4095. );
	hampSub->GetXaxis()->SetRangeUser(  -50.,   1000. );
	hint->GetXaxis()->SetRangeUser(     0.,  30000. );
	hintSub->GetXaxis()->SetRangeUser( -500.,  5000. );

//	hped->GetXaxis()->SetRangeUser(     0.,   500. );
//	hamp->GetXaxis()->SetRangeUser(     0.,   2000. );
//	hampSub->GetXaxis()->SetRangeUser(  -50.,   500. );
//	hint->GetXaxis()->SetRangeUser(     0.,  15000. );
//	hintSub->GetXaxis()->SetRangeUser( -500.,  2000. );

	  
	
	c1->Print( Form("data/fevme1_%06d.pdf(", vme_run), "pdf" );
	c2->Print( Form("data/fevme1_%06d.pdf", vme_run),  "pdf" );
	c3->Print( Form("data/fevme1_%06d.pdf", vme_run),  "pdf" );
	c4->Print( Form("data/fevme1_%06d.pdf", vme_run),  "pdf" );
	c5->Print( Form("data/fevme1_%06d.pdf)", vme_run), "pdf" );
	
	
	
//	fIn->Close();
	
	
	
	return;
}
