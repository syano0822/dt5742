// draw_event_2x8_overlay.C
//
// 2x8 view + all-channel overlay
// ┌─────────────────────┬─────────────────────────
// │ Canvas1 → Ch0–7 / Ch15–8 (2行×8列)            │
// │ Canvas2 → 16ch 全重ね (TLegend付き)           │
// └─────────────────────┴─────────────────────────
//
// 使い方：
//   root -l dt5742_float_ped.root
//   .L draw_event_2x8_overlay.C+
//   draw_event_2x8_overlay("dt5742_float_ped.root", 12, true);

#include "TFile.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include <iostream>

Color_t  colors[16] = {
    kBlack,        // 0
    kRed+1,        // 1
    kBlue+1,       // 2
    kGreen+2,      // 3
    kOrange+7,     // 4
    kMagenta+1,    // 5
    kCyan+1,       // 6
    kYellow+2,     // 7
    kViolet+1,     // 8
    kAzure+4,      // 9
    kPink+6,       // 10
    kTeal+2,       // 11
    kSpring+5,     // 12
    kGray+1,       // 13
    kRed-7,        // 14（ワイン色）
    kBlue-7        // 15（群青）
};

float min_disp_adc = 1200.;
float max_disp_adc = 4000.;

float min_disp_time = 150;
float max_disp_time = 175;

void draw_event_waveforms(const char* filename="wave_form.root",
			  int event=0, bool usePed=true)
{

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetStatBorderSize(0);
  gStyle->SetLegendBorderSize(0);

  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.07);

  
  TFile* fin=TFile::Open(filename,"READ");
    if(!fin||fin->IsZombie()){ std::cerr<<"File open error\n"; return; }

    char dname[64];
    sprintf(dname,"event_%06d",event);
    TDirectory* evtDir=(TDirectory*)fin->Get(dname);
    if(!evtDir){ std::cerr<<"No directory "<<dname<<"\n"; return; }

    TString tag = usePed? "ped":"raw";
    
    TCanvas* c1 = new TCanvas(Form("evt%d_view",event),
                              Form("Event %d (%s) 2x8 view",event,tag.Data()),
                              3200,1400);
    c1->Divide(8,2,0.002,0.005);

    for(int ch=0; ch<8; ch++){
        c1->cd(ch+1);
        TGraph* g=(TGraph*)evtDir->Get(Form("ch%02d_%s",ch,tag.Data()));
	g->SetLineColor(kRed+1);	
	g->SetLineWidth(2);
	g->SetMarkerStyle(20);
	g->SetMarkerSize(0.4);
	g->SetMaximum(0.02);
	g->SetMinimum(-0.2);
	g->SetMinimum(min_disp_adc);	
	g->SetMaximum(max_disp_adc);
	g->GetXaxis()->SetRangeUser(min_disp_time,max_disp_time);
        g->SetTitle(Form("Ch%02d",ch));
	gPad->SetGrid();
	g->GetXaxis()->SetTitleSize(0.06);
	g->GetYaxis()->SetTitleSize(0.06);
	g->GetXaxis()->SetTitleOffset(0.5);
	g->GetYaxis()->SetTitleOffset(1);
	g->Draw("ALP");
	
	TPaveText *pt = (TPaveText*)gPad->GetPrimitive("title");
	if (pt) {
	  pt->SetTextSize(0.05);
	}
	gPad->Modified();
	gPad->Update();
	
    }
    
    int pad=8;
    for(int ch=8; ch<=15; ch++){
        c1->cd(++pad);
        TGraph* g=(TGraph*)evtDir->Get(Form("ch%02d_%s",ch,tag.Data()));
	g->SetLineColor(kRed+1);	
	g->SetLineColor(kRed+1);	
	g->SetLineWidth(2);
	g->SetMarkerStyle(20);
	g->SetMarkerSize(0.4);
	g->SetMinimum(min_disp_adc);	
	g->SetMaximum(max_disp_adc);
	g->GetXaxis()->SetRangeUser(min_disp_time,max_disp_time);
        g->SetTitle(Form("Ch%02d",ch));
	gPad->SetGrid();
	g->GetXaxis()->SetTitleSize(0.06);
	g->GetYaxis()->SetTitleSize(0.06);
	g->GetXaxis()->SetTitleOffset(0.5);
	g->GetYaxis()->SetTitleOffset(1);
	g->Draw("ALP");
	
	TPaveText *pt = (TPaveText*)gPad->GetPrimitive("title");
	if (pt) {
	  pt->SetTextSize(0.05);
	}
	gPad->Modified();
	gPad->Update();

    }
    /*
      TCanvas* c2 = new TCanvas(Form("evt%d_overlay",event),
                              Form("Event %d (%s) overlay",event,tag.Data()),
                              1200,700);

    TLegend* leg = new TLegend(0.80,0.2,0.90,0.88);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);

    bool first=true;
    for(int ch=0; ch<16; ch++){
        TGraph* g=(TGraph*)evtDir->Get(Form("ch%02d_%s",ch,tag.Data()));
        if(!g) continue;
	
        Color_t c = 1+ch;               // 色分け
        g->SetLineColor(colors[ch]);
	g->SetLineWidth(2);
	g->SetMarkerStyle(20);
	g->SetMarkerSize(0.4);
	g->SetMinimum(min_disp_adc);	
	g->SetMaximum(max_disp_adc);
	g->GetXaxis()->SetRangeUser(min_disp_time,max_disp_time);
        g->SetTitle(Form("Ch%02d",ch));
	gPad->SetGrid();
	g->GetXaxis()->SetTitleSize(0.06);
	g->GetYaxis()->SetTitleSize(0.06);
	g->GetXaxis()->SetTitleOffset(0.5);
	g->GetYaxis()->SetTitleOffset(1);
        if(first){ g->Draw("AL"); first=false; }
        else       g->Draw("L SAME");
        leg->AddEntry(g,Form("Ch%02d",ch),"l");	
    }   
    leg->Draw();
    c2->Update();
    */
    
}
