///===================================================
/// PLOT EFFICIENCY
///===================================================
/// This method asumes the "reconstructed" histogram has the montecarlo bins
/// matched so that the division is alway lees than 1.
TGraphErrors* PlotEff(TH1F* montecarlo, 
                      TH1F* reconstructed, 
                      const char* title    ="", 
	                    const char* xtitle   ="", 
	                    const char* ytitle   ="", 
	                    const char* name     ="", 
	                    bool        logy     = 0,
	                    float       yMin     = 0, 
	                    float       yMax     = 0,
	                    const Char_t*     formula  = 0,
	                    float       r1       = 0, 
	                    float       r2       = 0, 
	                    TF1*        overlaid = 0,
	                    char*       label    = 0)
{
  //Int_t nbins1 = montecarlo ->GetNbinsX();
  //Int_t nbins2 = reconstructed->GetNbinsX();
  Int_t nbins1 = montecarlo ->GetXaxis()->GetLast();
  Int_t nbins2 = reconstructed->GetXaxis()->GetLast();
  if(nbins1!=nbins2) {
    cout << "ERROR, both histograms must have the same number of bins" << endl;
    cout << endl;
    return 0;
  } 
  const Int_t n = nbins1; 
  
  montecarlo->SetStats(kFALSE);
  reconstructed->SetStats(kFALSE);

  Float_t *x = new Float_t[n];
  Float_t *y = new Float_t[n];
  Float_t *ex = new Float_t[n];
  Float_t *ey = new Float_t[n];

  Float_t v1,v2;
  Float_t xBin = montecarlo->GetBinLowEdge(1) + 
    (montecarlo->GetBinWidth(1))/2;
  for(int i=1; i<=nbins1; i++) {
    v1 = montecarlo->GetBinContent(i);
    v2 = reconstructed->GetBinContent(i);
    if(v1!=0) {
      x[i-1] = xBin;
      y[i-1]  = v2/float(v1);
      ex[i-1] = montecarlo->GetBinWidth(i)/2;
      ey[i-1] = sqrt(y[i-1]*(1-y[i-1])/v1);
    } else {
      x[i-1] = xBin;
      y[i-1] = 0;
      ex[i-1] = 0;
      ey[i-1] = 0;
    }
    xBin += montecarlo->GetBinWidth(i);
  }
  //define graph
  TGraphErrors * result = new TGraphErrors(n,x,y,ex,ey);
  result->SetMarkerStyle(20);
  result->SetMarkerColor(kBlue);
  //result->SetMarkerSize(0.5);

  //define frame
  float xd = x[n-1]-x[0];
  float xmin = x[0] - 0.1*xd;
  float xmax = x[n-1] + 0.1*xd;

  //calculate y limit.
  // if it is not provided in the constructor, take 20 percent above
  // the maximum value.
  float ymax = 0;
  float ymin = 1000;
  if(yMax) {
    ymax = yMax;
  } else {
    for(int u=0; u<n;u++) {
      if(y[u]>ymax) ymax = y[u];
    }
    for(int u=0; u<n;u++) {
      if(y[u]<ymin) ymin = y[u];
    }
    float yd = ymax-ymin;
    ymin -= 0.1*yd;
    ymax += 0.2*yd;
  }
  

  TH2F* frame = new TH2F("frame","",2,xmin,xmax,2,yMin,ymax);
  frame->SetXTitle(xtitle);
  frame->SetYTitle(ytitle);
  frame->SetTitle(title);
  frame->SetStats(kFALSE);
  
  //plot
  TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
  if(!c1) {
    style();
    c1 = (TCanvas*)gROOT->FindObject("c1");
  }
  c1->SetLogy(0);
  c1->Clear();
  c1->SetBorderSize(2);
  if(logy) c1->SetLogy(1);
 
    
  //draw frame
  frame->Draw();
  c1->Modified();

  if(label) {
    TLatex* tex = new TLatex(0.2,0.85,label); tex->SetNDC();
    tex->SetTextSize(0.04); tex->SetTextColor(2); tex->SetLineWidth(2); 
    tex->Draw();
  }

  //draw eff.
  result->Draw("EP");
  if (formula) {
    TF1 *f1 = new TF1("f1",formula,xmin,xmax);
    if(r1==0 && r2==0) 
      result->Fit(f1);
    else
      result->Fit(f1,"","",r1,r2);
  }
  if(overlaid) overlaid->Draw("same");

  if(TString(name)!="")
    c1->SaveAs(name);
    
  delete frame;
  
  return result;
}