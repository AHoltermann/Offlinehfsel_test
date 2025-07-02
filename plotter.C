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

void plotter(){
  cout << "Start" << endl;
  
  cout << "Done!" << endl;
}
