/*
 * @Author: Zhiyong Lu(zhiyong.lu@cern.ch) 
 * @Date: 2022-10-03 14:46:45 
 * @Last Modified by: Zhiyong Lu
 * @Last Modified time: 2025-03-05 22:02:47
 */

#ifndef GRAPHSETTINGS
#define GRAPHSETTINGS
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "Rtypes.h"
#include "TPad.h"

const int color1s[] =
  {
    TColor::GetColor("#f94144"),
    TColor::GetColor("#f3722c"),
    TColor::GetColor("#f8961e"),
    TColor::GetColor("#f9844a"),
    TColor::GetColor("#f9c74f"),
    TColor::GetColor("#90be6d"),
    TColor::GetColor("#43aa8b"),
    TColor::GetColor("#4d908e"),
    TColor::GetColor("#577590"),
    TColor::GetColor("#277da1")
  };

const int color2s[] = {
  TColor::GetColor("#000000"),
  TColor::GetColor("#f26924"),
  TColor::GetColor("#9ccc3d"),
  TColor::GetColor("#15bdd4"),
  TColor::GetColor("#1a76bd"),
  TColor::GetColor("#6f4599"),
  TColor::GetColor("#ffd500"),
  TColor::GetColor("#bf1d89")
};

const int color3s[]={
    TColor::GetColor("#000000"),
    TColor::GetColor("#C11F7B"),
    TColor::GetColor("#CC0300"),
    TColor::GetColor("#EEBF13"),
    TColor::GetColor("#1D897C"),
    TColor::GetColor("#0055A3")
};

const int colorRed[]={
    TColor::GetColor("#dc2f02"),
    TColor::GetColor("#6a040f")
};

const int colorBlue[]={
    TColor::GetColor("#0096c7"),
    TColor::GetColor("#023e8a")
};

const int colorGreen[]={
    TColor::GetColor("#90a955"),
    TColor::GetColor("#31572c")
};

void SetMarkerAndLine(TGraphErrors* Graph, Int_t Color=0, Int_t style=0, Int_t linestyle=0, Double_t size=1){
    if(Color){
        Graph->SetMarkerColor(Color);
        Graph->SetLineColor(Color);
    }
    if(style){
        Graph->SetMarkerStyle(style);
        if(size){
            Graph->SetMarkerSize(size);
        }
    }
    if(linestyle){
        Graph->SetLineStyle(linestyle);
    }
}

void SetMarkerAndLine(TGraph* Graph, Int_t Color=0, Int_t style=0, Int_t linestyle=0, Double_t size=1){
    if(Color){
        Graph->SetMarkerColor(Color);
        Graph->SetLineColor(Color);
    }
    if(style){
        Graph->SetMarkerStyle(style);
        if(size){
            Graph->SetMarkerSize(size);
        }
    }
    if(linestyle){
        Graph->SetLineStyle(linestyle);
    }
}

void SetMarkerAndLine(TGraphAsymmErrors* Graph, Int_t Color=0, Int_t style=0, Int_t linestyle=0, Double_t size=1){
    if(Color){
        Graph->SetMarkerColor(Color);
        Graph->SetLineColor(Color);
    }
    if(style){
        Graph->SetMarkerStyle(style);
        if(size){
            Graph->SetMarkerSize(size);
        }
    }
    if(linestyle){
        Graph->SetLineStyle(linestyle);
    }
}

void SetMarkerAndLine(TH1D* Graph, Int_t Color=0, Int_t style=0, Int_t linestyle=0, Double_t size=1){
    if(Color){
        Graph->SetMarkerColor(Color);
        Graph->SetLineColor(Color);
    }
    if(style){
        Graph->SetMarkerStyle(style);
        if(size){
            Graph->SetMarkerSize(size);
        }
    }
    if(linestyle){
        Graph->SetLineStyle(linestyle);
    }
}

void DrawTwoPad(TCanvas* canvas, Double_t y, std::vector<TGraphErrors*>& UpGraph, TLegend* UpLegend, std::vector<Int_t>& UpColor, std::vector<Int_t>& UpMarkerStyle, std::vector<Int_t>& UpLineStyle, std::vector<TGraphErrors*>& DownGraph, TLegend* DownLegend, std::vector<Int_t>& DownColor, std::vector<Int_t>& DownMarkerStyle, std::vector<Int_t>& DownLineStyle, Double_t Upxmin=0, Double_t Upxmax=0, Double_t Upymin=0, Double_t Upymax=0, Double_t Downxmin=0, Double_t Downxmax=0, Double_t Downymin=0, Double_t Downymax=0){
    //UpPad
    canvas->cd();
    TPad* pad1 = new TPad("pad1", "pad1", 0, y, 1, 1);
    pad1->Draw();
	pad1->cd();
	pad1->SetBorderMode(0);
	pad1->SetBorderSize(2);
	pad1->SetLeftMargin(0.11);
	pad1->SetRightMargin(0.01);
	pad1->SetTopMargin(0.1);
	pad1->SetBottomMargin(0);
	pad1->SetTicks(1,1);
    for(Int_t i=0;i<UpGraph.size();i++){
        UpGraph[i]->SetMarkerColor(UpColor[i]);
        UpGraph[i]->SetLineColor(UpColor[i]);
        UpGraph[i]->SetMarkerStyle(UpMarkerStyle[i]);
        UpGraph[i]->SetLineStyle(UpLineStyle[i]);
    }
    if(!Upxmin&&!Upxmax)UpGraph[0]->GetXaxis()->SetRangeUser(Upxmin,Upxmax);
    if(!Upymin&&!Upymax)UpGraph[0]->GetYaxis()->SetRangeUser(Upymin,Upymax);
    UpGraph[0]->Draw("AP");
    for(Int_t i=1;i<UpGraph.size();i++)UpGraph[i]->Draw("P");
    UpLegend->Draw();


    //DownPad
    canvas->cd();
    TPad* pad2 = new TPad("pad2", "pad2", 0, 0, 1, y);
    pad2->Draw();
	pad2->cd();
	pad2->SetBorderMode(0);
	pad2->SetBorderSize(2);
	pad2->SetLeftMargin(0.08);
	pad2->SetRightMargin(0.01);
	pad2->SetTopMargin(0.0);
	pad2->SetBottomMargin(0.2);
	pad2->SetTicks(1,1);
    for(Int_t i=0;i<DownGraph.size();i++){
        DownGraph[i]->SetMarkerColor(DownColor[i]);
        DownGraph[i]->SetLineColor(DownColor[i]);
        DownGraph[i]->SetMarkerStyle(DownMarkerStyle[i]);
        DownGraph[i]->SetLineStyle(DownLineStyle[i]);
    }
    if(!Downxmin&&!Downxmax)DownGraph[0]->GetXaxis()->SetRangeUser(Downxmin,Downxmax);
    if(!Downymin&&!Downymax)DownGraph[0]->GetYaxis()->SetRangeUser(Downymin,Downymax);
    DownGraph[0]->Draw("AP");
    for(Int_t i=1;i<DownGraph.size();i++)DownGraph[i]->Draw("P");
    DownLegend->Draw();
}


void ShiftAlongXaxis(TGraphErrors *ge, Double_t shift)
{
    // Shift original TGraphErrors along x-axis by Double_t shift.
    
    if(!ge){std::cout<<"!ge"<<std::endl;exit(0);}
    
    Int_t nPoints = ge->GetN();
    Double_t x = 0.;
    Double_t y = 0.;
    for(Int_t p=0;p<nPoints;p++)
    {
        ge->GetPoint(p,x,y);
        x+=shift;
        ge->SetPoint(p,x,y);
    } // end of for(Int_t p=0;p<nPoints;p++)
} // end of void ShiftAlongXaxis(TGraphErrors *ge, Double_t shift)



void myGraphSetUp(TGraphErrors *currentGraph=0, Float_t currentMarkerSize = 1.0,
                  int currentMarkerStyle=20, int currentMarkerColor=0,
                  int currentLineStyle=1, int currentLineColor=0, int currentFillColor=0, int currentFillStyle=0)
{
    currentGraph->SetMarkerSize(currentMarkerSize);
    currentGraph->SetMarkerStyle(currentMarkerStyle);
    currentGraph->SetMarkerColor(currentMarkerColor);
    currentGraph->SetLineStyle(currentLineStyle);
    currentGraph->SetLineColor(currentLineColor);
    currentGraph->SetFillColor(currentFillColor);
    currentGraph->SetFillStyle(currentFillStyle);
    return;
}

void myPadSetUp(TPad *currentPad, float currentLeft, float currentTop, float currentRight, float currentBottom){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
}


void SetStyle(Bool_t graypalette=false) {
  std::cout << "Setting style!" << std::endl;
  
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);


}

void FakeHistosOnlyForExample(TH1* &hstat, TH1* &hsyst, TH1*&hsystCorr) {

  TF1 * fExpo = new TF1 ("fExpo", "expo");
  fExpo->SetParameters(10, -0.3);
  hstat     = new TH1F("hstat", "hstat", 100, 0, 10);
  hsyst     = new TH1F("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1F("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  hstat->FillRandom("fExpo",20000);
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
			  hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
			  hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
			  hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }
  

}

TCanvas* newCanvas(char* name, char* title, Int_t ww, Int_t wh){
    TCanvas *c = new TCanvas(name, title, ww, wh);
    c->Range(0,0,1,1);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    // c->SetLeftMargin(0.03333334);
    // c->SetRightMargin(0.01);
    // c->SetTopMargin(0.04761905);
    // c->SetBottomMargin(0.05);
    c->SetFrameBorderMode(0);
    return c;
}

void SetGraphLabel(TGraph* vn){
    vn->GetXaxis()->SetLabelSize(0.08);
    vn->GetXaxis()->SetTitleSize(0.08);
    vn->GetXaxis()->SetTitleOffset(0.8);
    vn->GetYaxis()->SetLabelSize(0.05);
    vn->GetYaxis()->SetTitleSize(0.07);
    vn->GetYaxis()->SetTitleOffset(0.9);
}

void SetCanvas(TCanvas* canvas1){
    canvas1->Range(0,0,1,1);
    canvas1->SetFillColor(0);
    canvas1->SetBorderMode(0);
    canvas1->SetBorderSize(2);
    canvas1->SetLeftMargin(0.03333334);
    canvas1->SetRightMargin(0.01);
    canvas1->SetTopMargin(0.04761905);
    canvas1->SetBottomMargin(0.05);
    canvas1->SetFrameBorderMode(0);
}

void DrawPad(TPad* pad, Double_t LeftMargin, Double_t BottomMargin, Double_t TopMargin=0.){
    pad->Draw();
    pad->cd();
    pad->Range(0.,0.06, 6.,0.4);
    pad->SetFillColor(0);
    pad->SetBorderMode(0);
    pad->SetBorderSize(2);
    pad->SetLeftMargin(LeftMargin);
    pad->SetRightMargin(0.05);
    pad->SetTopMargin(TopMargin);
    pad->SetBottomMargin(BottomMargin);
    pad->SetFrameBorderMode(0);
}

void DrawOne(double xmin, double xmax){
    TGraph* One = new TGraph();
    One->SetPoint(0,xmin,1);
    One->SetPoint(1,xmax,1);
    One->SetLineColor(kBlack);
    One->SetLineStyle(kDashed);
    One->Draw("L");
}

#endif
