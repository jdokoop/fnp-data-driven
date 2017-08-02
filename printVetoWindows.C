#include <iostream>

using namespace std;

void printVetoWindows()
{
	// parameters for pT dependent veto cut
	//    f(pT) = p[0] + p[1]/(p[2]+1)
	static const double par[4][3] = {
		{0.005, 0.05, -0.3}, // B0
		{0.005, 0.05, -0.3}, // B1
		{0.010, 0.08, -0.3}, // B2
		{0.018, 0.08, -0.3}  // B3
	};

	// fixed limits for veto cut taken from the original values.
	// pT dependent cut is too wide for low pT to kill most of the tracks
	// which is not desired.
	// For high pT, the pT dep cut is too small since the peak at 0 can be made
	// by mis-association to the cluster from the pairs. This mis-association
	// is made in the negative side of the peak. therefore, the veto cut used
	// in the negative side is just taken.
	static const double veto_limit[4][2] = {
		{0.04, 0.02},  // B0
		{0.06, 0.02},  // B1
		{0.08, 0.04},  // B2
		{0.08, 0.02},  // B3
	};

	static double ptcut[4][2] = {{ -1, -1}, { -1, -1}, { -1, -1}, { -1, -1}};

	static bool init_flag = false;
	if (!init_flag) {
		for (int ilay = 0; ilay < 4; ilay++)
		{
			for (int i = 0; i < 2; i++)
			{
				ptcut[ilay][i] = par[ilay][1] / (veto_limit[ilay][i] - par[ilay][0]) - par[ilay][2];
				cout << "pT cut : " << ilay << " " << i << " " << ptcut[ilay][i] << endl;
			}
		}
	}

	cout << endl;
	cout << "For ptcut0 < pT < ptcut1: " << endl;
	for (int layer = 0; layer < 4; layer++)
	{
		cout << "   B" << layer << ": " <<  par[layer][0] <<  " + " << par[layer][1] << " / ( pT + " << par[layer][2] << ")"  << endl;
	}

	cout << endl;
	double limit;
	double pt = 1.5;
	int layer = 1;
	if (pt < ptcut[layer][0])      { limit = veto_limit[layer][0]; /*cout<<"lim_h "<<pt<<endl;*/}
	else if (ptcut[layer][1] < pt) { limit = veto_limit[layer][1]; /*cout<<"lim_l "<<pt<<endl;*/}
	else                        { limit = par[layer][0] + par[layer][1] / (pt + par[layer][2]); /*cout<<"lim_pt "<<pt<<endl;*/}
	cout << "LIMIT (" << pt << ", " << layer << ") = " << limit << endl << endl;

	//Plot example window for B0
	TF1 *f_win[4];

	for(int layer=0; layer<4; layer++)
	{
		f_win[layer] = new TF1(Form("f_win_%i",layer), Form("(x < %g)*0.04 + (x >= %g && x < %g)*(%g + %g / ( x + %g)) + (x >= %g)*0.02", ptcut[layer][0], ptcut[layer][0], ptcut[layer][1], par[layer][0], par[layer][1], par[layer][2], ptcut[layer][1]), 0.0, 9.0);
		cout << Form("(x < %g)*0.04 + (x >= %g && x < %g)*(%g + %g / ( x + %g)) + (x >= %g)*0.02", ptcut[layer][0], ptcut[layer][0], ptcut[layer][1], par[layer][0], par[layer][1], par[layer][2], ptcut[layer][1]) << endl;
	}

	TCanvas *cB0 = new TCanvas("cB0", "cB0", 600, 600);
	f_win[0]->SetTitle("");
	f_win[0]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	f_win[0]->GetXaxis()->SetRangeUser(0.0, 9.0);
	f_win[0]->GetYaxis()->SetTitle("charge #times #Delta#phi Upper Limit");
	f_win[0]->GetYaxis()->SetTitleOffset(1.4);
	f_win[0]->GetYaxis()->SetRangeUser(0.0, 0.08);
	f_win[0]->Draw();
	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.07);
	latex.DrawLatex(.15, .8, "B0");

	TCanvas *cB1 = new TCanvas("cB1", "cB1", 600, 600);
	f_win[1]->SetTitle("");
	f_win[1]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	f_win[1]->GetXaxis()->SetRangeUser(0.0, 9.0);
	f_win[1]->GetYaxis()->SetTitle("charge #times #Delta#phi Upper Limit");
	f_win[1]->GetYaxis()->SetTitleOffset(1.4);
	f_win[1]->GetYaxis()->SetRangeUser(0.0, 0.08);
	f_win[1]->Draw();
	latex.DrawLatex(.15, .8, "B1");

	TCanvas *cB2 = new TCanvas("cB2", "cB2", 600, 600);
	f_win[2]->SetTitle("");
	f_win[2]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	f_win[2]->GetXaxis()->SetRangeUser(0.0, 9.0);
	f_win[2]->GetYaxis()->SetTitle("charge #times #Delta#phi Upper Limit");
	f_win[2]->GetYaxis()->SetTitleOffset(1.4);
	f_win[2]->GetYaxis()->SetRangeUser(0.0, 0.08);
	f_win[2]->Draw();
	latex.DrawLatex(.15, .8, "B2");

	TCanvas *cB3 = new TCanvas("cB3", "cB3", 600, 600);
	f_win[3]->SetTitle("");
	f_win[3]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	f_win[3]->GetXaxis()->SetRangeUser(0.0, 9.0);
	f_win[3]->GetYaxis()->SetTitle("charge #times #Delta#phi Upper Limit");
	f_win[3]->GetYaxis()->SetTitleOffset(1.4);
	f_win[3]->GetYaxis()->SetRangeUser(0.0, 0.08);
	f_win[3]->Draw();
	latex.DrawLatex(.15, .8, "B3");
}