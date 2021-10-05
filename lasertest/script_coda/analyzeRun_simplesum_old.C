void analyzeRun_simplesum(int vme_run) {
	
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(111111);	
	
	const int n=6;

	int nped_a[n] = {2,2,2,2,2,2};
//	int nped_a[n] = {12,12,12,12,12,12};
//	int nped_a[n] = {6,6,12,12,12,12};
	int nped_b[n] = {20,20,20,20,20,20};
	int nint_a[n] = {20,20,20,20,20,20};
//	int nint_b[n] = {35,35,35,35,35,35};
//	int nint_b[n] = {48,48,48,48,48,48};
//	int nint_b[n] = {47,47,47,47,47,47};
	int nint_b[n] = {95,95,95,95,95,95};

	int nped[n];
	
	// peak around bin , full width about  bin
	int nint[n];

	for( int i = 0; i < n; i++ ) {
		nped[i] = nped_b[i] - nped_a[i];
		nint[i] = nint_b[i] - nint_a[i];
	}

	TH1I *hampSub_total= new TH1I( "hampSub_total", "Amplitude (pedestal subtracted)", 10500, -500., 10000. );
	TH1I *hintSub_total= new TH1I( "hintSub_total", "Integral (pedestal subtracted)", 30500, -500., 30000. );

	TH2I *hampSub_2D= new TH2I( "hampSub_2D", "Amplitude (pedestal subtracted);sum64;quadsum", 4000, -500., 3500., 4000, -500., 3500.);
	TH2I *hintSub_2D= new TH2I( "hintSub_2D", "Integral (pedestal subtracted);sum64;quadsum", 8500, -500., 8000., 8500, -500., 8000.);

	TH1I *hped[n],*hamp[n],*hampSub[n],*hint[n],*hintSub[n];

	for( int i = 0; i < n; i++ ) {
	hped[i] = new TH1I( Form("hped%i",i), "Pedestal",   4096, 0., 4096. );
	hped[i]->SetTitle( Form("Pedestal Distribution (average fADC count between %d - %d ns)", 
		nped_a[i]*4, nped_b[i]*4) );
	hped[i]->GetXaxis()->SetTitle( "Pedestal [fADC counts]" );
	
	hamp[i]    = new TH1I( Form("hamp%i",i),    "Amplitude",                       4096, 0., 4096. );
	hampSub[i] = new TH1I( Form("hampSub%i",i), "Amplitude (pedestal subtracted)", 5000, -1000., 4000. );
	hamp[i]->SetTitle( "Amplitude Distribution (no pedestal subtraction)" );
	hamp[i]->GetXaxis()->SetTitle( "Amplitude [fADC Counts]" );
	hampSub[i]->SetTitle( "Amplitude Distribution (with pedestal subtraction)" );
	hampSub[i]->GetXaxis()->SetTitle( "Amplitude [fADC]" );
	
	hint[i]    = new TH1I( Form("hint%i",i),    "Integral",                       101000, -1000., 100000. );
	hintSub[i] = new TH1I( Form("hintSub%i",i), "Integral (pedestal subtracted)", 31000, -1000., 30000. );
	hint[i]->SetTitle( "Integral Distribution (no pedestal subtraction)" );
	hint[i]->GetXaxis()->SetTitle( "Integral [fADC Counts]" );
	hintSub[i]->SetTitle( "Integral Distribution (with pedestal subtraction)" );
	hintSub[i]->GetXaxis()->SetTitle( "Integral [fADC]" );
	}

	TFile *fIn = new TFile( Form("data/fevme1_%06d.root",vme_run), "READ" );
	TTree *T   = (TTree*)fIn->Get( "T" );
	
	int nwave;
	int wave[10][100];
	T->SetBranchAddress( "wave_n", &nwave );
	T->SetBranchAddress( "wave"  , &wave  );
		
	cout << "N of events " << T->GetEntries() << endl;
	for( int ievt = 2; ievt < T->GetEntries()-11; ievt++ ) {
	  
	  
	  if( ievt%100000 == 0 ) cout << "Processing event " << ievt << endl;
	  T->GetEvent(ievt);
	  
	  int ampSub_total=0,intSub_total=0;
	  int ampSub_sum64=0,intSub_sum64=0;
	  for( int i = 0; i < n; i++ ) {
	    
	    //	if (i==1) continue;	
	    int loc_ped = 0;
	    if (i<2) for( int isamp = nped_a[i]; isamp < nped_b[i]; isamp++ ) loc_ped += wave[i+3][isamp];
	    else for( int isamp = nped_a[i]; isamp < nped_b[i]; isamp++ ) loc_ped += 4095-wave[i+3][isamp];
	    
	    double loc_amp = 0.; double loc_int = 0.;
	    for( int isamp = nint_a[i]; isamp < nint_b[i]; isamp++ ) {	    
	      double f_wave;
	      if (i<2) f_wave = static_cast<double>( wave[i+3][isamp] );
	      else  f_wave = static_cast<double>( 4095 - wave[i+3][isamp] );
	      // f_wave = static_cast<double>( wave[i+3][isamp] );
	      loc_int += f_wave;
	      if( f_wave > loc_amp ) loc_amp = f_wave;
	      // loc_amp = f_wave;
	      // hamp[i]->Fill( loc_amp );
	    }
	    
	    double f_ped = static_cast<double>( loc_ped ) / static_cast<double>( nped[i] );
	    
	    hped[i]->Fill( f_ped );
	    hamp[i]->Fill( loc_amp );
	    hint[i]->Fill( loc_int );
	    
	    
	    loc_amp = loc_amp - f_ped;
	    loc_int = loc_int - f_ped*static_cast<double>( nint[i] );
	    
	    hampSub[i]->Fill( loc_amp );
	    hintSub[i]->Fill( loc_int );
	    
	    if (2<=i && i<=5) {
	      ampSub_total += loc_amp;
	      intSub_total += loc_int;
	    }
	    
	    if (i==0) {ampSub_sum64=loc_amp; intSub_sum64=loc_int;}
	    
	  }
	  
	  hampSub_total->Fill( ampSub_total );
	  hintSub_total->Fill( intSub_total );
	  
	  hampSub_2D->Fill(ampSub_total,ampSub_sum64 );
	  hintSub_2D->Fill(intSub_total,intSub_sum64 );
	  
	} // end of event loop

	hampSub_total->SetLineColor( 7 );
	hintSub_total->SetLineColor( 7 );
	hampSub_total->SetLineStyle( 2 );
	hintSub_total->SetLineStyle( 2 );
	hampSub_total->SetLineWidth( 2 );
	hintSub_total->SetLineWidth( 2 );
	hampSub_total->Rebin( 4 );
	hintSub_total->Rebin( 4 );


	for( int i = 0; i < n; i++ ) {
	hamp[i]->Rebin(4);
	hampSub[i]->Rebin(4);
	hint[i]->Rebin(4);
	hintSub[i]->Rebin(4);
	
	hped[i]->SetLineWidth( 2 );
	hamp[i]->SetLineWidth( 2 );
	hampSub[i]->SetLineWidth( 2 );
	hint[i]->SetLineWidth( 2 );
	hintSub[i]->SetLineWidth( 2 );
	
	hped[i]->SetLineColor( i+1 );
	hamp[i]->SetLineColor( i+1 );
	hampSub[i]->SetLineColor( i+1 );
	hint[i]->SetLineColor( i+1 );
	hintSub[i]->SetLineColor( i+1 );

	}
	
	char *label[n]={"sum64","dynode","quad1","quad2","quad3","quad4"};
	TLegend *leg = new TLegend(0.75, 0.4-0.03*n, 0.95, 0.4);
	for(int i=0;i<n;i++){
  	leg->AddEntry(hped[i],label[i], "l");  
	}
  	leg->AddEntry(hampSub_total,"quadsum", "l");  

	double max=0; 
	
	
	TCanvas *c1 = new TCanvas( "c1", "ped", 1000, 1000 );
	c1->SetTickx(); c1->SetTicky();
	c1->cd(); 
	for( int i = 0; i < n; i++ ) hped[i]->Draw("same");
	for( int i = 0; i < n; i++ ) {if (hped[i]->GetMaximum() > max) max=hped[i]->GetMaximum();}
	hped[0]->SetMaximum(max);
	leg->Draw();
	
	TCanvas *c2 = new TCanvas( "c2", "amp", 1000, 1000 );
	c2->SetTickx(); c2->SetTicky();
	c2->cd(); 
	for( int i = 0; i < n; i++ ) hamp[i]->Draw("same");
	max=0; 
	for( int i = 0; i < n; i++ ) {if (hamp[i]->GetMaximum() > max) max=hamp[i]->GetMaximum();}
	hamp[0]->SetMaximum(max);
	leg->Draw();
	
	TCanvas *c3 = new TCanvas( "c3", "int", 1000, 1000 );
	c3->SetTickx(); c3->SetTicky();
	c3->cd(); 
	for( int i = 0; i < n; i++ ) hint[i]->Draw("same");
	max=0; 
	for( int i = 0; i < n; i++ ) {if (hint[i]->GetMaximum() > max) max=hint[i]->GetMaximum();}
	hint[0]->SetMaximum(max);
	leg->Draw();
	
	int index=0; //fit sum64
//	int index=1; //fit dynode
	 double cmax, xmax;
	 double lowend=300,width=50; 

	TCanvas *c4 = new TCanvas( "c4", "ampSub", 1000, 1000 );
	c4->SetTickx(); c4->SetTicky();
	c4->cd(); 
	for( int i = 0; i < n; i++ ) hampSub[i]->Draw("same");
	max=0; 
	for( int i = 0; i < n; i++ ) {if (i==1) continue; if (hampSub[i]->GetMaximum() > max) max=hampSub[i]->GetMaximum();}
	hampSub[0]->SetMaximum(max);
	leg->Draw();
	hampSub_total->Draw("same");
	gPad->SetLogy(1);
	
	  TF1 *f1 = new TF1( "f1", "gaus", 0., 1. );
	  hampSub[index]->GetXaxis()->SetRangeUser( lowend, 4000. );
	  cmax = static_cast<double>( hampSub[index]->GetMaximum() );
	  xmax = static_cast<double>( hampSub[index]->GetBinCenter(hampSub[index]->GetMaximumBin()) );
	  f1->SetRange( xmax-(xmax-width), xmax+(xmax-width) );
	  f1->SetParameters( cmax, xmax, (xmax-width)/2. );
	  hampSub[index]->Fit( "f1", "R0Q" );
	  width = f1->GetParameter(2);
	  f1->SetRange( xmax-width*1.5, xmax+width*1.5 );
	  hampSub[index]->Fit( "f1", "R0", "same" );
	
	width = 50;
	TCanvas *c5 = new TCanvas( "c5", "intSub", 1000, 1000 );
	c5->SetTickx(); c5->SetTicky();
	c5->cd(); 
	for( int i = 0; i < n; i++ ) hintSub[i]->Draw("same");
	max=0; 
	for( int i = 0; i < n; i++ ) {if (i==1) continue; if (hintSub[i]->GetMaximum() > max) max=hintSub[i]->GetMaximum();}
	hintSub[0]->SetMaximum(max);
	leg->Draw();
	hintSub_total->Draw("same");
	gPad->SetLogy(1);
	
	  TF1 *f2 = new TF1( "f2", "gaus", 0., 1. );
	  hintSub[index]->GetXaxis()->SetRangeUser( lowend, 10000. );
	  cmax = static_cast<double>( hintSub[index]->GetMaximum() );
	  xmax = static_cast<double>( hintSub[index]->GetBinCenter(hintSub[index]->GetMaximumBin()) );
	  f2->SetRange( xmax-(xmax-width*2), xmax+(xmax-width*2) );
	  f2->SetParameters( cmax, xmax, (xmax-width*2)/2. );	  
	  hintSub[index]->Fit( "f2", "R0Q" );
	  width = f2->GetParameter(2);
	  f2->SetRange( xmax-width*1.5, xmax+width*1.5 );
	  hintSub[index]->Fit( "f2", "R0", "same" );
	
	TCanvas *c6 = new TCanvas( "c6", "amp_2D", 1000, 1000 );
	c6->SetTickx(); c6->SetTicky();
	c6->cd();
	hampSub_2D->Draw("colz");

	TCanvas *c7 = new TCanvas( "c7", "int_2D", 1000, 1000 );
	c7->SetTickx(); c7->SetTicky();
	c7->cd();
	hintSub_2D->Draw("colz");
	
	
	for( int i = 0; i < n; i++ ) {
	hped[i]->GetXaxis()->SetRangeUser(     1000,   3000. );
	hamp[i]->GetXaxis()->SetRangeUser(     0.,   4000. );
	hampSub[i]->GetXaxis()->SetRangeUser( -200.,   5000. );
	hint[i]->GetXaxis()->SetRangeUser(    40000.,  80000. );
	hintSub[i]->GetXaxis()->SetRangeUser( -200.,  10000. );
	}
	hampSub_total->GetXaxis()->SetRangeUser( -200.,  1000. );
	hintSub_total->GetXaxis()->SetRangeUser( -200.,  2500. );	
	
	hintSub_2D->GetXaxis()->SetRangeUser(-200.,8000.);
	hintSub_2D->GetYaxis()->SetRangeUser(-200.,8000.);
	
	
	c1->Print( Form("data/fevme1_%06d.pdf(", vme_run), "pdf" );
	c2->Print( Form("data/fevme1_%06d.pdf", vme_run),  "pdf" );
	c3->Print( Form("data/fevme1_%06d.pdf", vme_run),  "pdf" );
	c4->Print( Form("data/fevme1_%06d.pdf", vme_run),  "pdf" );
	c5->Print( Form("data/fevme1_%06d.pdf", vme_run), "pdf" );
	c6->Print( Form("data/fevme1_%06d.pdf", vme_run), "pdf" );
	c7->Print( Form("data/fevme1_%06d.pdf)", vme_run), "pdf" );
	
	fIn->Close();
		
	return;
}
