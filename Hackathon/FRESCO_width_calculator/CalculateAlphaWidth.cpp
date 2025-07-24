#include <TGraph.h>
#include <TCanvas.h>

double Sa = 10.61474; //MeV
double mu_26Mg = 22.*4./26.;

int A1 = 4, A2 = 22, Z1 = 2, Z2 = 10;

double prefactor = pow(TMath::Hbar() * TMath::C() / TMath::Qe() / 1.e6 * 1.e15,2.) / (2 * mu_26Mg * 938);

double Ex = 10.949; //MeV

int L = 1;

int main()
{
    TGraph *graph = new TGraph();
    graph->SetName("PartialWidth");
    int counter = 0;

    ifstream input;
    input.open("fort.46");

    char line[256];
    for(int j=0;j<1;j++) input.getline(line, sizeof(line));

    double radius, u_r, a, b, c, d, wf2;

    double Ecm = (Ex - Sa); //I may need to adjust because I can't remember if my code takes CoM energy or lab energy...

    while(input >> radius >> a >> u_r >> b >> c>> d)
    {
        wf2 = pow(u_r/radius,2.);
        double reducedWidth = prefactor * radius * wf2 * 1e3;
        double partialWidth = 2 * reducedWidth * penetrability(Z1, Z2, A1, A2, 1, Ecm, radius);
//         double penetrability(int Z1, int Z2, double A1, double A2, int L, double E, double r)
        graph->SetPoint(counter,radius,partialWidth*1000);
        counter++;
    }

    graph->SaveAs("PartialWidthGraph.root");
    TCanvas *c1 = new TCanvas("c1","c1",800,800);
    graph->Draw("ALP");
    graph->GetYaxis()->SetTitle("#Gamma_{sp} / eV"); graph->GetYaxis()->CenterTitle(); graph->GetYaxis()->SetTitleOffset(1.1);
    graph->GetXaxis()->SetTitle("R / fm"); graph->GetXaxis()->CenterTitle();
    c1->SetLogy();
    c1->SaveAs("PartialWidthGraph.png");
    c1->SaveAs("PartialWidthGraph.C");
}
