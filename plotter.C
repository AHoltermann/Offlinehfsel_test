using namespace std;

void Plotter1D(vector<TObject*> objects,  // these vectors should all be the same size, with any TF1 objects plotted last
               vector<int> colors, 
               vector<string> labels,
               const char* xAxisTitle,
               const char* yAxisTitle,
               float xMin = 0.0,
               float xMax = 10.0,
               float yMin = 0.0,
               float yMax = 100.0,
               float legendX = 0.7,
               float legendY = 0.7,
               int logy = 0,
               const char* outputFileName = "plotter_output.pdf") {


    TCanvas *c = new TCanvas("c", "Canvas", 800, 600);
    TLegend* legend = new TLegend(legendX, legendY, legendX + 0.18, legendY + 0.18);
    int colorIdx = 0;
    int labelIdx = 0;
    for (size_t i = 0; i < objects.size(); ++i) {
        TH1D *h = dynamic_cast<TH1D*>(objects[i]);
        TF1  *f = dynamic_cast<TF1*>(objects[i]);
        if (h) {
            h->SetLineColor(colors[colorIdx]);
            h->SetTitle(labels[labelIdx].c_str());
            h->GetXaxis()->SetTitle(xAxisTitle);
            h->GetYaxis()->SetTitle(yAxisTitle);
            h->GetXaxis()->SetRangeUser(xMin, xMax);
            h->GetYaxis()->SetRangeUser(yMin, yMax);
            if (logy) c->SetLogy();
            h->Draw(colorIdx == 0 ? "" : "SAME");
            legend->AddEntry(h, labels[labelIdx].c_str(), "l");
            colorIdx++;
            labelIdx++;
        } else if (f) {
            f->SetLineColor(colors[colorIdx]);
            f->SetTitle(labels[labelIdx].c_str());
            if (logy) c->SetLogy();
            f->Draw(colorIdx == 0 ? "" : "SAME");
            legend->AddEntry(f, labels[labelIdx].c_str(), "l");
            colorIdx++;
            labelIdx++;
        }
    }
    legend->Draw();

    c->SaveAs(outputFileName);

}


//function to draw percent lines on ratio plot
void DrawPercents(TF1 *func, double percent, float range_low, float range_high, Color_t color, TLegend *legend) 
{
    TF1 *f_diff0 = new TF1("f_diff0", Form("func(x) - %f", percent), range_low, range_high);
    double x_cross = f_diff0->GetX(0, range_low, range_high);

    TLine *vline = new TLine(x_cross, 0, x_cross, percent);
    vline->SetLineColor(color);
    vline->SetLineStyle(2);
    vline->Draw("SAME");

    TLine *hline = new TLine(0, percent, x_cross, percent);
    hline->SetLineColor(color);
    hline->SetLineStyle(2);
    hline->Draw("SAME");

    legend->AddEntry(vline, Form("%.0f%% at x = %.2f", percent*100, x_cross), "l");
}

void RatioPlot(TH1D *h,
               const char* title,
               const char* xAxisTitle ,
               const char* yAxisTitle, 
               const char* outputFileName = "plotter_output.pdf", 
               int ptCut = 0,
               //vector<float> init_guess,,
               float xMin = 0.0,
               float xMax = 100.0,
               float yMin = 0.0,
               float yMax = 1.01,
               float legendX = 0.65,
               float legendY = 0.15,
               float fit_low = 5,
               float fit_high = 200.0,                
               vector<int> colors = {kBlack, kRed, kBlue, kGreen+2, kRed},
               vector<float> percents = {0.9, 0.95, 0.99}) {


    TCanvas *c = new TCanvas("c", "Canvas", 800, 600);
    TLegend* legend = new TLegend(legendX, legendY, legendX + 0.25, legendY + 0.25);
    legend->SetTextSize(0.03);
    int colorIdx = 0;
    int labelIdx = 0;
    h->SetLineColor(colors[colorIdx]);
    h->SetTitle(title);
    h->GetXaxis()->SetTitle(xAxisTitle);
    h->GetYaxis()->SetTitle(yAxisTitle);
    h->GetXaxis()->SetRangeUser(xMin, xMax);
    h->GetYaxis()->SetRangeUser(yMin, yMax);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(1);
    h->Draw("P");
    legend->AddEntry(h, "HIJING OO", "p");
    colorIdx++;
    labelIdx++;
    
    TF1 *func = new TF1("func", "[0]*exp([1]*x) + [2]*exp([3]*x) + [4]", fit_low, fit_high);
    func->SetParameters(-9.98, -0.498, 73.32, -1.30, 1.00);
    //func->SetParameters(-0.5, -0.0633, -4.64, -0.469, 0.99);
    h->Fit("func", "0", "", fit_low, fit_high);

    func->SetLineColor(colors[colorIdx]);
    func->SetLineWidth(2);
    func->SetLineStyle(1);
    func->Draw("SAME");
    legend->AddEntry(func, "Fit", "l");
    colorIdx++;
    labelIdx++;

    DrawPercents(func, percents[0], fit_low, fit_high, colors[colorIdx], legend);
    DrawPercents(func, percents[1], fit_low, fit_high, colors[colorIdx + 1], legend);
    DrawPercents(func, percents[2], fit_low, fit_high, colors[colorIdx + 2], legend);
    
    legend->Draw();

    TLatex *tex = new TLatex();
    tex->SetTextSize(0.03);
    tex->SetNDC();
    tex->DrawLatex(0.40, 0.35, "PVFilter");
    tex->DrawLatex(0.40, 0.30, "CCFilter");
    tex->DrawLatex(0.40, 0.25, "|VZ| < 15");
    tex->DrawLatex(0.40, 0.20, "nVtx > 0");
    if (ptCut != 0)
        tex->DrawLatex(0.40, 0.15, Form("track p_{T} > %i GeV/c", ptCut));

    c->SaveAs(outputFileName);

}

void plotter(){
    cout << "Start" << endl;
    gStyle->SetOptStat(0);
    TFile* f = new TFile("online.root");

    RatioPlot((TH1D*)f->Get("ratio14AND"), "Ratio of events with and without online 14 AND cut",
    "minimum HFEMax [GeV]", "Ratio", "ratio14AND.pdf");
    RatioPlot((TH1D*)f->Get("ratio14ANDpt2"),"Ratio of events with and without online 14 AND, pt > 2 cut",
    "minimum HFEMax [GeV]", "Ratio", "ratio14ANDpt2.pdf", 2);
    RatioPlot((TH1D*)f->Get("ratio14ANDpt3"),"Ratio of events with and without online 14 AND, pt > 3 cut",
    "minimum HFEMax [GeV]", "Ratio", "ratio14ANDpt3.pdf", 3);
    RatioPlot((TH1D*)f->Get("ratio14ANDpt4"),"Ratio of events with and without online 14 AND, pt > 4 cut",
    "minimum HFEMax [GeV]", "Ratio", "ratio14ANDpt4.pdf", 4);

    RatioPlot((TH1D*)f->Get("ratio14OR"),"Ratio of events with and without online 14 OR cut",
    "maximum HFEMax [GeV]", "Ratio", "ratio14OR.pdf");
    RatioPlot((TH1D*)f->Get("ratio14ORpt2"),"Ratio of events with and without online 14 OR, pt > 2 cut",
    "maximum HFEMax [GeV]", "Ratio", "ratio14ORpt2.pdf", 2);
    RatioPlot((TH1D*)f->Get("ratio14ORpt3"),"Ratio of events with and without online 14 OR, pt > 3 cut",
    "maximum HFEMax [GeV]", "Ratio", "ratio14ORpt3.pdf", 3);
    RatioPlot((TH1D*)f->Get("ratio14ORpt4"),"Ratio of events with and without online 14 OR, pt > 4 cut",
    "maximum HFEMax [GeV]", "Ratio", "ratio14ORpt4.pdf", 4);

    RatioPlot((TH1D*)f->Get("ratio16AND"),"Ratio of events with and without online 16 AND cut",
    "minimum HFEMax [GeV]", "Ratio", "ratio16AND.pdf");
    RatioPlot((TH1D*)f->Get("ratio16ANDpt2"),"Ratio of events with and without online 16 AND, pt > 2 cut",
    "minimum HFEMax [GeV]", "Ratio", "ratio16ANDpt2.pdf", 2);
    RatioPlot((TH1D*)f->Get("ratio16ANDpt3"),"Ratio of events with and without online 16 AND, pt > 3 cut",
    "minimum HFEMax [GeV]", "Ratio", "ratio16ANDpt3.pdf", 3);
    RatioPlot((TH1D*)f->Get("ratio16ANDpt4"),"Ratio of events with and without online 16 AND, pt > 4 cut",
    "minimum HFEMax [GeV]", "Ratio", "ratio16ANDpt4.pdf", 4);

    RatioPlot((TH1D*)f->Get("ratio16OR"),"Ratio of events with and without online 16 OR cut",
    "maximum HFEMax [GeV]", "Ratio", "ratio16OR.pdf");
    RatioPlot((TH1D*)f->Get("ratio16ORpt2"),"Ratio of events with and without online 16 OR, pt > 2 cut",
    "maximum HFEMax [GeV]", "Ratio", "ratio16ORpt2.pdf", 2);
    RatioPlot((TH1D*)f->Get("ratio16ORpt3"),"Ratio of events with and without online 16 OR, pt > 3 cut",
    "maximum HFEMax [GeV]", "Ratio", "ratio16ORpt3.pdf", 3);
    RatioPlot((TH1D*)f->Get("ratio16ORpt4"),"Ratio of events with and without online 16 OR, pt > 4 cut",
    "maximum HFEMax [GeV]", "Ratio", "ratio16ORpt4.pdf", 4);
  cout << "Done!" << endl;
}
