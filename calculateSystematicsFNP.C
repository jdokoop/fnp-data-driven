#include <iostream>

using namespace std;

//----------------------------------
// Variables
//----------------------------------

//Save plots?
bool savePlots = false;

//Compute errors using TEfficiency?
bool effErrors = true;

//Ratio of etas and pizeros to photons
float etaPhotonRatio = 1.0;
float pizeroPhotonRatio = 1.0;

//Which particle variations to group curves by?
//0: Photon
//1: Pion
//2: Eta
int colorCodeVariation = 0;

//Rebin factor
const int REBINF = 40;

const int NBINSFNP = 9;

//Number of variations for each particle species
const int NVARETAS = 5;
const int NVARPIZEROS = 3;
const int NVARPHOTONS = 4;
const int NVARETAPIRATIO = 3;

//Eta pi ratio variations
float etaPiRatio[NVARETAPIRATIO] = {0.96, 1.04, 0.88};

//Spectra
TF1 *f_photon_spectrum[NVARPHOTONS];
TF1 *f_pizero_spectrum[NVARPIZEROS];
TF1 *f_eta_spectrum[NVARETAS];

//Fit to FNP with veto cut applied
TF1 *f_fnp_fit_with_veto;

//Fit to FNP without veto cut applied
TF1 *f_fnp_fit_without_veto;

//Total evaluated FNPs
const int NUMFNPS = NVARPIZEROS * NVARPHOTONS * NVARETAS * NVARETAPIRATIO;

//RMS of all systematic variations on FNP
float rmsFNP_noveto[NBINSFNP];
float rmsFNP_noveto_asymm_above[NBINSFNP];
float rmsFNP_noveto_asymm_below[NBINSFNP];

float rmsFNP_veto[NBINSFNP];
float rmsFNP_veto_asymm_above[NBINSFNP];
float rmsFNP_veto_asymm_below[NBINSFNP];

TH1D *h_rms_plus;
TH1D *h_rms_minus;

//TGraph to draw systematics as polygon
TGraph *g_systematicsRMS_noveto;
TGraph *g_systematicsRMS_veto;

//Simulated electrons
TH1D *h_elec_pT_pizeros[NVARPIZEROS];
TH1D *h_elec_pT_pizeros_rejected[NVARPIZEROS];
TH1D *h_elec_pT_pizeros_accepted[NVARPIZEROS];

TH1D *h_elec_pT_photons[NVARPHOTONS];
TH1D *h_elec_pT_photons_rejected[NVARPHOTONS];
TH1D *h_elec_pT_photons_accepted[NVARPHOTONS];

TH1D *h_elec_pT_etas[NVARETAS];
TH1D *h_elec_pT_etas_rejected[NVARETAS];
TH1D *h_elec_pT_etas_accepted[NVARETAS];

//Data electrons + hadrons
TH1D *h_elec_pT;
TH1D *h_elec_pT_rejected;
TH1D *h_elec_pT_accepted;

TH1D *h_elec_sw_pT;
TH1D *h_elec_sw_pT_rejected;
TH1D *h_elec_sw_pT_accepted;

TH1D *h_hadron_pT;
TH1D *h_hadron_pT_rejected;
TH1D *h_hadron_pT_accepted;

//Inclusive electrons (contamination subtracted)
TH1D *h_elec_pT_inclusive;
TH1D *h_elec_pT_isolated;

//Survival rate for each particle species
TH1D *h_survival_photons[NVARPHOTONS];
TH1D *h_survival_pizeros[NVARPIZEROS];
TH1D *h_survival_etas[NVARETAS];

//Ingredients for FNP calculation
TH1D *h_conversion_efficiency[NUMFNPS];
TH1D *h_killed_efficiency;

//FNP
TH1D *h_P[NUMFNPS];
TH1D *h_NP[NUMFNPS];
TH1D *h_FNP[NUMFNPS];
TGraphAsymmErrors *g_eff_FNP[NUMFNPS];

//FNP with conversion veto applied
TH1D *h_P_conv[NUMFNPS];
TH1D *h_NP_conv[NUMFNPS];
TH1D *h_FNP_conv[NUMFNPS];
TGraphAsymmErrors *g_eff_FNP_conv[NUMFNPS];

//Histos to draw the ratio of all systematic variations to the default FNP
TH1D *h_FNP_default_ratios[NUMFNPS];
TH1D *h_FNP_conv_default_ratios[NUMFNPS];

//TBoxes for systematics
TBox *b_fnp_syst_noveto[NBINSFNP];
TBox *b_fnp_syst_veto[NBINSFNP];

//----------------------------------
// Functions
//----------------------------------

void formatHistograms(TH1D *& h, string xTitle, string yTitle, string title)
{
	h->SetTitle(title.c_str());

	h->GetXaxis()->SetTitleFont(62);
	h->GetXaxis()->SetLabelFont(62);
	h->GetXaxis()->SetTitleOffset(1.3);
	h->GetXaxis()->SetTitle(xTitle.c_str());

	h->GetYaxis()->SetTitleFont(62);
	h->GetYaxis()->SetLabelFont(62);
	h->GetYaxis()->SetTitleOffset(1.3);
	h->GetYaxis()->SetTitle(yTitle.c_str());

	h->SetMarkerStyle(20);
	h->SetMarkerSize(0.6);
	h->SetMarkerColor(kAzure + 2);
	h->SetLineColor(kAzure + 2);
}

void formatGraphs(TGraphAsymmErrors *& h, string xTitle, string yTitle, string title)
{
	h->SetTitle(title.c_str());

	h->GetXaxis()->SetTitleFont(62);
	h->GetXaxis()->SetLabelFont(62);
	h->GetXaxis()->SetTitleOffset(1.3);
	h->GetXaxis()->SetTitle(xTitle.c_str());

	h->GetYaxis()->SetTitleFont(62);
	h->GetYaxis()->SetLabelFont(62);
	h->GetYaxis()->SetTitleOffset(1.3);
	h->GetYaxis()->SetTitle(yTitle.c_str());

	h->SetMarkerStyle(20);
	h->SetMarkerSize(0.6);
	h->SetMarkerColor(kAzure + 2);
	h->SetLineColor(kAzure + 2);
}


void readFiles()
{
	//Files from simulation
	TFile *fSimsPizero = new TFile("Sims/Cocktail062617_2/twopizeros_cdphiana.root");
	TFile *fSimsPhoton = new TFile("Sims/Cocktail062617_2/twophotons_cdphiana.root");
	TFile *fSimsEta = new TFile("Sims/Cocktail062617_2/twoetas_cdphiana.root");

	for (int i = 0; i < NVARPIZEROS; i++)
	{
		h_elec_pT_pizeros[i]          = (TH1D*) fSimsPizero->Get(Form("h_scaled_elec_pT_%i", i));
		h_elec_pT_pizeros_rejected[i] = (TH1D*) fSimsPizero->Get(Form("h_scaled_elec_pT_rejected_pTcut_%i", i));
		h_elec_pT_pizeros_accepted[i] = (TH1D*) fSimsPizero->Get(Form("h_scaled_elec_pT_notrejected_pTcut_%i", i));
	}

	for (int i = 0; i < NVARETAS; i++)
	{
		h_elec_pT_etas[i]          = (TH1D*) fSimsEta->Get(Form("h_scaled_elec_pT_%i", i));
		h_elec_pT_etas_rejected[i] = (TH1D*) fSimsEta->Get(Form("h_scaled_elec_pT_rejected_pTcut_%i", i));
		h_elec_pT_etas_accepted[i] = (TH1D*) fSimsEta->Get(Form("h_scaled_elec_pT_notrejected_pTcut_%i", i));
	}

	for (int i = 0; i < NVARPHOTONS; i++)
	{
		h_elec_pT_photons[i]          = (TH1D*) fSimsPhoton->Get(Form("h_scaled_elec_pT_%i", i));
		h_elec_pT_photons_rejected[i] = (TH1D*) fSimsPhoton->Get(Form("h_scaled_elec_pT_rejected_pTcut_%i", i));
		h_elec_pT_photons_accepted[i] = (TH1D*) fSimsPhoton->Get(Form("h_scaled_elec_pT_notrejected_pTcut_%i", i));
	}

	//Files from data
	TFile *fData = new TFile("Data/11167/calcsurvivalrate.root");

	h_elec_pT          = (TH1D*) fData->Get("h_elec_pT");
	h_elec_pT_rejected = (TH1D*) fData->Get("h_elec_pT_rejected_pTcut");
	h_elec_pT_accepted = (TH1D*) fData->Get("h_elec_pT_notrejected_pTcut");

	h_elec_sw_pT          = (TH1D*) fData->Get("h_elec_sw_pT");
	h_elec_sw_pT_rejected = (TH1D*) fData->Get("h_elec_sw_pT_rejected_pTcut");
	h_elec_sw_pT_accepted = (TH1D*) fData->Get("h_elec_sw_pT_notrejected_pTcut");

	h_hadron_pT          = (TH1D*) fData->Get("h_hadron_pT");
	h_hadron_pT_rejected = (TH1D*) fData->Get("h_hadron_pT_rejected_pTcut");
	h_hadron_pT_accepted = (TH1D*) fData->Get("h_hadron_pT_notrejected_pTcut");
}


void setSimErrorsToZero()
{
	//Set the errors on histograms from simulation to zero.
	//This is done to determine the extent to which the statistical errors on the final FNP
	//are driven by having low statistics in simulations (which can be easily fixed!) or in data...

	for (int i = 0; i < NVARPIZEROS; i++)
	{
		for (int j = 1; j <= h_elec_pT_pizeros[i]->GetNbinsX(); j++)
		{
			h_elec_pT_pizeros[i]->SetBinError(j, 0.0);
			h_elec_pT_pizeros_rejected[i]->SetBinError(j, 0.0);
			h_elec_pT_pizeros_accepted[i]->SetBinError(j, 0.0);
		}
	}

	for (int i = 0; i < NVARETAS; i++)
	{
		for (int j = 1; j <= h_elec_pT_etas[i]->GetNbinsX(); j++)
		{
			h_elec_pT_etas[i]->SetBinError(j, 0.0);
			h_elec_pT_etas_rejected[i]->SetBinError(j, 0.0);
			h_elec_pT_etas_accepted[i]->SetBinError(j, 0.0);
		}
	}

	for (int i = 0; i < NVARPHOTONS; i++)
	{
		for (int j = 1; j <= h_elec_pT_photons[i]->GetNbinsX(); j++)
		{
			h_elec_pT_photons[i]->SetBinError(j, 0.0);
			h_elec_pT_photons_rejected[i]->SetBinError(j, 0.0);
			h_elec_pT_photons_accepted[i]->SetBinError(j, 0.0);
		}
	}
}


float getEtaPhotonRatio(int etaIndex, int photonIndex)
{
	//Take the ratio of the total yield of etas to photons for a given set of spectra
	float etaYield = f_eta_spectrum[etaIndex]->Eval(18.0);
	float photonYield = f_photon_spectrum[photonIndex]->Eval(18.0);

	return etaYield / photonYield;
}


float getPizeroPhotonRatio(int pizeroIndex, int photonIndex)
{
	//Take the ratio of the total yield of etas to photons for a given set of spectra
	float pizeroYield = f_pizero_spectrum[pizeroIndex]->Eval(18.0);
	float photonYield = f_photon_spectrum[photonIndex]->Eval(18.0);

	return pizeroYield / photonYield;
}


void rebinHistograms()
{
	/*
	h_elec_pT_pizeros->Rebin(REBINF);
	h_elec_pT_pizeros_rejected->Rebin(REBINF);
	h_elec_pT_pizeros_accepted->Rebin(REBINF);
	h_elec_pT_photons->Rebin(REBINF);
	h_elec_pT_photons_rejected->Rebin(REBINF);
	h_elec_pT_photons_accepted->Rebin(REBINF);
	h_elec_pT_etas->Rebin(REBINF);
	h_elec_pT_etas_rejected->Rebin(REBINF);
	h_elec_pT_etas_accepted->Rebin(REBINF);

	h_elec_pT->Rebin(REBINF);
	h_elec_pT_rejected->Rebin(REBINF);
	h_elec_pT_accepted->Rebin(REBINF);
	h_elec_sw_pT->Rebin(REBINF);
	h_elec_sw_pT_rejected->Rebin(REBINF);
	h_elec_sw_pT_accepted->Rebin(REBINF);
	h_hadron_pT->Rebin(REBINF);
	h_hadron_pT_rejected->Rebin(REBINF);
	h_hadron_pT_accepted->Rebin(REBINF);
	*/


	//Rebin data electron spectrum with variable bin width
	double bins[10] = {1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};

	for (int i = 0; i < NVARPIZEROS; i++)
	{
		h_elec_pT_pizeros[i] = (TH1D*) h_elec_pT_pizeros[i]->Rebin(NBINSFNP, Form("h_elec_pT_pizeros_rebinned_%i", i), bins);
		h_elec_pT_pizeros_rejected[i] = (TH1D*) h_elec_pT_pizeros_rejected[i]->Rebin(NBINSFNP, Form("h_elec_pT_pizeros_rejected_rebinned_%i", i), bins);
		h_elec_pT_pizeros_accepted[i] = (TH1D*) h_elec_pT_pizeros_accepted[i]->Rebin(NBINSFNP, Form("h_elec_pT_pizeros_accepted_rebinned_%i", i), bins);
	}

	for (int i = 0; i < NVARETAS; i++)
	{
		h_elec_pT_etas[i] = (TH1D*) h_elec_pT_etas[i]->Rebin(NBINSFNP, Form("h_elec_pT_etas_rebinned_%i", i), bins);
		h_elec_pT_etas_rejected[i] = (TH1D*) h_elec_pT_etas_rejected[i]->Rebin(NBINSFNP, Form("h_elec_pT_etas_rejected_rebinned_%i", i), bins);
		h_elec_pT_etas_accepted[i] = (TH1D*) h_elec_pT_etas_accepted[i]->Rebin(NBINSFNP, Form("h_elec_pT_etas_accepted_rebinned_%i", i), bins);
	}

	for (int i = 0; i < NVARPHOTONS; i++)
	{
		h_elec_pT_photons[i] = (TH1D*) h_elec_pT_photons[i]->Rebin(NBINSFNP, Form("h_elec_pT_photons_rebinned_%i", i), bins);
		h_elec_pT_photons_rejected[i] = (TH1D*) h_elec_pT_photons_rejected[i]->Rebin(NBINSFNP, Form("h_elec_pT_photons_rejected_rebinned_%i", i), bins);
		h_elec_pT_photons_accepted[i] = (TH1D*) h_elec_pT_photons_accepted[i]->Rebin(NBINSFNP, Form("h_elec_pT_photons_accepted_rebinned_%i", i), bins);
	}

	h_elec_pT = (TH1D*) h_elec_pT->Rebin(NBINSFNP, "h_elec_pT_rebinned", bins);
	h_elec_pT_rejected = (TH1D*) h_elec_pT_rejected->Rebin(NBINSFNP, "h_elec_pT_rejected_rebinned", bins);
	h_elec_pT_accepted = (TH1D*) h_elec_pT_accepted->Rebin(NBINSFNP, "h_elec_pT_accepted_rebinned", bins);
	h_elec_sw_pT = (TH1D*) h_elec_sw_pT->Rebin(NBINSFNP, "h_elec_sw_pT_rebinned", bins);
	h_elec_sw_pT_rejected = (TH1D*) h_elec_sw_pT_rejected->Rebin(NBINSFNP, "h_elec_sw_pT_rejected_rebinned", bins);
	h_elec_sw_pT_accepted = (TH1D*) h_elec_sw_pT_accepted->Rebin(NBINSFNP, "h_elec_sw_pT_accepted_rebinned", bins);
	h_hadron_pT = (TH1D*) h_hadron_pT->Rebin(NBINSFNP, "h_hadron_pT_rebinned", bins);
	h_hadron_pT_rejected = (TH1D*) h_hadron_pT_rejected->Rebin(NBINSFNP, "h_hadron_pT_rejected_rebinned", bins);
	h_hadron_pT_accepted = (TH1D*) h_hadron_pT_accepted->Rebin(NBINSFNP, "h_hadron_pT_accepted_rebinned", bins);
}


void defineSpectra()
{
	f_photon_spectrum[0] = new TF1("f_photon_spectrum_0", "2*TMath::Pi()*x*[0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])", 0.0, 18.0);
	f_photon_spectrum[0]->SetParameters(0.242345, -0.0827585, 0.00918447, 4.13943, 13.6974);

	f_photon_spectrum[1] = new TF1("f_photon_spectrum_1", "2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
	f_photon_spectrum[1]->SetParameters(13.0488, -0.191588, 0.0164036, 0.999159, 8.42105);

	f_photon_spectrum[2] = new TF1("f_photon_spectrum_2", "2*TMath::Pi()*x*(([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])))", 0.0, 18.0);
	f_photon_spectrum[2]->SetParameters(0.171216, -0.0312528, 0.00442746, 11.6472, 29.538);

	f_photon_spectrum[3] = new TF1("f_photon_spectrum_3", "2*TMath::Pi()*x*(([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])))", 0.0, 18.0);
	f_photon_spectrum[3]->SetParameters(0.012215, -0.0539112, 0.00527051, 4.69066, 11.7111);


	f_pizero_spectrum[0] = new TF1("f_pizero_spectrum_0", "2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
	f_pizero_spectrum[0]->SetParameters(254.066, 0.470186, 0.0380066, 0.737713, 8.28442);

	f_pizero_spectrum[1] = new TF1("f_pizero_spectrum_1", "2*TMath::Pi()*x*(([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])))", 0.0, 18.0);
	f_pizero_spectrum[1]->SetParameters(339.852, 0.417738, 0.0599019, 0.719301, 8.33941);

	f_pizero_spectrum[2] = new TF1("f_pizero_spectrum_2", "2*TMath::Pi()*x*(([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])))", 0.0, 18.0);
	f_pizero_spectrum[2]->SetParameters(184.594, 0.524929, 0.0111792, 0.762548, 8.24103);


	f_eta_spectrum[0] = new TF1("f_eta_spectrum_0", "TMath::Power(TMath::Sqrt(1+(0.135/x)*(0.135/x)),-1.0)*0.5*TMath::Sqrt(1+(0.135/x)*(0.135/x))*(2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])))", 0.0, 18.0);
	f_eta_spectrum[0]->SetParameters(254.066, 0.470186, 0.0380066, 0.737713, 8.28442);

	f_eta_spectrum[1] = new TF1("f_eta_spectrum_1", "2*TMath::Pi()*x*(([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])))", 0.0, 18.0);
	f_eta_spectrum[1]->SetParameters(1201.68, -0.0135258, 0.0215599, 0.782018, 9.31206);

	f_eta_spectrum[2] = new TF1("f_eta_spectrum_2", "2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
	f_eta_spectrum[2]->SetParameters(1372.01, 0.0103744, 0.0225786, 0.71981, 9.07938);

	f_eta_spectrum[3] = new TF1("f_eta_spectrum_3", "2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
	f_eta_spectrum[3]->SetParameters(545.206, 0.421156, -0.00156003, 0.58791, 8.04994);

	f_eta_spectrum[4] = new TF1("f_eta_spectrum_4", "2*TMath::Pi()*x*(([0]*(([1]-1)*([1]-1))/(([1]*[2]+0.53*([1]-1))*([1]*[2]+0.53))*pow(([1]*[2]+TMath::Sqrt(0.53*0.53+x*x))/([1]*[2]+0.53),-1*[1])))", 0.0, 18.0);
	f_eta_spectrum[4]->SetParameters(0.786558, 9.24328, 0.100848);
}


void calculateSpeciesSurvival()
{
	//Survival = Accepted / Total
	for (int i = 0; i < NVARETAS; i++)
	{
		h_survival_etas[i] = (TH1D*) h_elec_pT_etas_accepted[i]->Clone(Form("h_survival_etas_%i", i));
		h_survival_etas[i]->Divide(h_elec_pT_etas[i]);
	}

	for (int i = 0; i < NVARPIZEROS; i++)
	{
		h_survival_pizeros[i] = (TH1D*) h_elec_pT_pizeros_accepted[i]->Clone(Form("h_survival_pizeros_%i", i));
		h_survival_pizeros[i]->Divide(h_elec_pT_pizeros[i]);
	}

	for (int i = 0; i < NVARPHOTONS; i++)
	{
		h_survival_photons[i] = (TH1D*) h_elec_pT_photons_accepted[i]->Clone(Form("h_survival_photons_%i", i));
		h_survival_photons[i]->Divide(h_elec_pT_photons[i]);
	}
}


void calculateConversionEfficiency()
{
	//Conversion efficiency = photonic electrons after cut / photonic electrons without cut
	//Here we calculate the conversion efficiency for all possible combinations of reweighting each particle spectrum
	//There are 60 total combinations, labeled by three indices (i, j, k) = (photon, pizero, eta)

	//cout << "********* Indices = (photon, pizero, eta, eta-pi ratio) **********" << endl << endl;

	TH1D *h_pT_photonic[NUMFNPS];
	TH1D *h_pT_photonic_accepted[NUMFNPS];

	int numHisto = 0;
	for (int i = 0; i < NVARPHOTONS; i++)
	{
		for (int j = 0; j < NVARPIZEROS; j++)
		{
			for (int k = 0; k < NVARETAS; k++)
			{
				for (int m = 0; m < NVARETAPIRATIO; m++)
				{
					//Scale photonic electron yields to account for relative ratio between pions (etas) and photons
					float pionScale = getPizeroPhotonRatio(j, i);
					float etaScale = getEtaPhotonRatio(i, k);

					//Scale electrons from eta by variation on eta/pizero ratio
					float etaPiRatioScale = etaPiRatio[m];

					TH1D *h_aux_pizeros = (TH1D*) h_elec_pT_pizeros[j]->Clone(Form("h_aux_pizeros_%i_%i_%i_%i", i, j, k, m));
					TH1D *h_aux_etas = (TH1D*) h_elec_pT_etas[k]->Clone(Form("h_aux_etas_%i_%i_%i_%i", i, j, k, m));

					h_aux_pizeros->Scale(pionScale);
					h_aux_etas->Scale(etaScale);
					h_aux_etas->Scale(etaPiRatioScale);

					h_pT_photonic[numHisto] = (TH1D*) h_elec_pT_photons[i]->Clone(Form("h_pT_photonic_%i_%i_%i_%i", i, j, k, m));
					h_pT_photonic[numHisto]->Add(h_aux_pizeros);
					h_pT_photonic[numHisto]->Add(h_aux_etas);

					//cout << numHisto << ": " << i << " " << j << " " << k << " " << m << endl;

					numHisto++;
				}
			}
		}
	}

	numHisto = 0;
	for (int i = 0; i < NVARPHOTONS; i++)
	{
		for (int j = 0; j < NVARPIZEROS; j++)
		{
			for (int k = 0; k < NVARETAS; k++)
			{
				for (int m = 0; m < NVARETAPIRATIO; m++)
				{
					//Scale photonic electron yields to account for relative ratio between pions (etas) and photons
					float pionScale = getPizeroPhotonRatio(j, i);
					float etaScale = getEtaPhotonRatio(i, k);

					//Scale electrons from eta by variation on eta/pizero ratio
					float etaPiRatioScale = etaPiRatio[m];

					TH1D *h_aux_pizeros = (TH1D*) h_elec_pT_pizeros_accepted[j]->Clone(Form("h_aux_pizeros_accepted_%i_%i_%i_%i", i, j, k, m));
					TH1D *h_aux_etas = (TH1D*) h_elec_pT_etas_accepted[k]->Clone(Form("h_aux_etas_accepted_%i_%i_%i_%i", i, j, k, m));

					h_aux_pizeros->Scale(pionScale);
					h_aux_etas->Scale(etaScale);
					h_aux_etas->Scale(etaPiRatioScale);

					h_pT_photonic_accepted[numHisto] = (TH1D*) h_elec_pT_photons_accepted[i]->Clone(Form("h_pT_photonic_%i_%i_%i_%i", i, j, k, m));
					h_pT_photonic_accepted[numHisto]->Add(h_aux_pizeros);
					h_pT_photonic_accepted[numHisto]->Add(h_aux_etas);
					numHisto++;
				}
			}
		}
	}

	for (int i = 0; i < NUMFNPS; i++)
	{
		h_conversion_efficiency[i] = (TH1D*) h_pT_photonic_accepted[i]->Clone(Form("h_conversion_efficiency_%i", i));
		h_conversion_efficiency[i]->Divide(h_pT_photonic[i]);
	}
}


void calculateRandomlyKilledEfficiency()
{
	//Randomly killed efficiency = hadrons after cut / hadrons without cut
	h_killed_efficiency = (TH1D*) h_hadron_pT_accepted->Clone("h_killed_efficiency");
	h_killed_efficiency->Divide(h_hadron_pT);
}


void getCleanElectronSample()
{
	//Subtract swapped tracks from electron candidates
	h_elec_pT_inclusive = (TH1D*) h_elec_pT->Clone("h_elec_pT_inclusive");
	h_elec_pT_inclusive->Add(h_elec_sw_pT, -1.0);

	h_elec_pT_isolated = (TH1D*) h_elec_pT_accepted->Clone("h_elec_pT_isolated");
	h_elec_pT_isolated->Add(h_elec_sw_pT_accepted, -1.0);
}


int getPhotonIndex(int index)
{
	return index / 45;
}


int getPionIndex(int index)
{
	return (index - (index / 45) * 45) / 15;
}


int getEtaIndex(int index)
{
	return (index - (index / 15) * 15) / 3;
}

int getEtaPiRatioIndex(int index)
{
	return index % 3;
}


void getFNP()
{
	for (int i = 0; i < NUMFNPS; i++)
	{
		//Histogram filled with all ones
		TH1D *h_one = (TH1D*) h_elec_pT_inclusive->Clone("h_one");
		for (int i = 1; i <= h_one->GetNbinsX(); i++)
		{
			h_one->SetBinContent(i, 1.0);
			h_one->SetBinError(i, 0.0);
		}

		//Calculate P and NP without conv veto
		TH1D *h_aux1 = (TH1D*) h_elec_pT_inclusive->Clone("h_aux1");
		h_aux1->Multiply(h_killed_efficiency);

		TH1D *h_aux2 = (TH1D*) h_elec_pT_isolated->Clone("h_aux2");
		h_aux2->Add(h_aux1, -1.0);

		TH1D *h_aux3 = (TH1D*) h_conversion_efficiency[i]->Clone("h_aux3");
		h_aux3->Add(h_one, -1.0);

		TH1D *h_aux4 = (TH1D*) h_killed_efficiency->Clone("h_aux4");
		h_aux4->Multiply(h_aux3);

		h_P[i] = (TH1D*) h_aux2->Clone("h_P");
		h_P[i]->Divide(h_aux4);

		h_NP[i] = (TH1D*) h_elec_pT_inclusive->Clone("h_NP");
		h_NP[i]->Add(h_P[i], -1.0);

		//Calculate FNP without conv veto
		TH1D *h_aux5 = (TH1D*) h_NP[i]->Clone("h_aux5");
		h_aux5->Add(h_P[i]);

		h_FNP[i] = (TH1D*) h_NP[i]->Clone("h_FNP");
		h_FNP[i]->Divide(h_aux5);

		//Calculate FNP using BayesDivide to handle the errors,
		//so they don't go above 1
		g_eff_FNP[i] = new TGraphAsymmErrors();
		g_eff_FNP[i]->Divide(h_NP[i], h_aux5, "cl=0.683 b(1,1) mode");

		//Set horizontal error bar on TGraphAsymmErrors to zero
		for (int j = 0; j < g_eff_FNP[i]->GetN(); j++)
		{
			g_eff_FNP[i]->SetPointError(j, 0.0, 0.0, g_eff_FNP[i]->GetErrorYlow(j), g_eff_FNP[i]->GetErrorYhigh(j));
		}

		//Calculate FNP with the conversion verto by multiplying NP by the randomly killed
		//effiency, and P by the isolation cut efficiency and randomly killed efficiency
		//Then, FNP = e_r*NP / (e_r*NP + e_r*e_i*P)
		TH1D *h_aux6 = (TH1D*) h_NP[i]->Clone("h_aux6");
		TH1D *h_aux7 = (TH1D*) h_P[i]->Clone("h_aux7");

		h_aux6->Multiply(h_killed_efficiency);
		h_aux7->Multiply(h_killed_efficiency);
		h_aux7->Multiply(h_conversion_efficiency[i]);

		TH1D *h_aux8 = (TH1D*) h_aux6->Clone("h_aux8");
		h_aux8->Add(h_aux7);

		h_FNP_conv[i] = (TH1D*) h_aux6->Clone("h_FNP_conv");
		h_FNP_conv[i]->Divide(h_aux8);

		//Calculate FNP using BayesDivide class to handle the errors,
		//so they don't go above 1
		g_eff_FNP_conv[i] = new TGraphAsymmErrors();
		g_eff_FNP_conv[i] = new TGraphAsymmErrors();
		g_eff_FNP_conv[i]->Divide(h_aux6, h_aux8, "cl=0.683 b(1,1) mode");

		//Set horizontal error bar on TGraphAsymmErrors to zero
		for (int j = 0; j < g_eff_FNP_conv[i]->GetN(); j++)
		{
			g_eff_FNP_conv[i]->SetPointError(j, 0.0, 0.0, g_eff_FNP_conv[i]->GetErrorYlow(j), g_eff_FNP_conv[i]->GetErrorYhigh(j));
		}
	}
}


void getAsymmRMS()
{
	//----------------------------------------------------------------
	// RMS Calculation for FNP Without Veto Cut
	//----------------------------------------------------------------

	//Display the systematics as a filled polygon
	float x_syst[NBINSFNP * 2];
	float y_syst[NBINSFNP * 2];

	h_rms_plus = (TH1D*) h_FNP[0]->Clone("h_rms_plus");
	h_rms_minus = (TH1D*) h_FNP[0]->Clone("h_rms_minus");

	//Get the RMS of all systematic variations at every pT point
	for (int i = 1; i <= NBINSFNP; i++)
	{
		//Vectors to store, for every FNP bin, the values of the systematic variations that are above and below the default value of FNP
		std::vector<float> valuesAbove;
		std::vector<float> valuesBelow;

		//Array to store, for every FNP bin, the ratio of the systematic variation to the default value
		float valuesRatio[NUMFNPS] = {0.0};
		float values[NUMFNPS] = {0.0};

		for (int j = 0; j < NUMFNPS; j++)
		{
			float fnpDefault = h_FNP[0]->GetBinContent(i);
			float fnp = h_FNP[j]->GetBinContent(i);

			if (fnp > fnpDefault)
			{
				valuesAbove.push_back(fnp - fnpDefault);
			}
			else if (fnp < fnpDefault)
			{
				valuesBelow.push_back(fnp - fnpDefault);
			}

			values[j] = fnp;
			valuesRatio[j] = fnp / h_FNP[0]->GetBinContent(i);
		}

		//Calculate RMS for variations above the default
		rmsFNP_noveto_asymm_above[i - 1] = 0.0;
		for (int k = 0; k < valuesAbove.size(); k++)
		{
			rmsFNP_noveto_asymm_above[i - 1] += valuesAbove[k] * valuesAbove[k];
		}
		rmsFNP_noveto_asymm_above[i - 1] = TMath::Sqrt((float) rmsFNP_noveto_asymm_above[i - 1] / valuesAbove.size());

		//Calculate RMS for variations below the default
		rmsFNP_noveto_asymm_below[i - 1] = 0.0;
		for (int k = 0; k < valuesBelow.size(); k++)
		{
			rmsFNP_noveto_asymm_below[i - 1] += valuesBelow[k] * valuesBelow[k];
		}
		rmsFNP_noveto_asymm_below[i - 1] = TMath::Sqrt((float) rmsFNP_noveto_asymm_below[i - 1] / valuesBelow.size());

		float mean = TMath::Mean(NUMFNPS, values);
		float meanRatio = TMath::Mean(NUMFNPS, valuesRatio);
		float rmsFNP_novetoRatio = TMath::RMS(NUMFNPS, valuesRatio);
		h_rms_plus->SetBinContent(i, meanRatio + rmsFNP_noveto_asymm_above[i - 1]/mean);
		h_rms_minus->SetBinContent(i, meanRatio - rmsFNP_noveto_asymm_below[i - 1]/mean);

		b_fnp_syst_noveto[i - 1] = new TBox(h_FNP[0]->GetBinCenter(i) - 0.15, h_FNP[0]->GetBinContent(i) - rmsFNP_noveto_asymm_below[i - 1], h_FNP[0]->GetBinCenter(i) + 0.15, h_FNP[0]->GetBinContent(i) + rmsFNP_noveto_asymm_above[i - 1]);
	}

	//Fill TGraph
	//For a given pT, there will be 2 entries in the TGraph:
	//Central value \pm RMS
	for (int i = 1; i <= 2 * NBINSFNP; i++)
	{
		x_syst[i - 1] = h_FNP[0]->GetBinCenter((i - 1) % NBINSFNP + 1);

		if (i - 1 < NBINSFNP)
		{
			y_syst[i - 1] = h_FNP[0]->GetBinContent((i - 1) % NBINSFNP + 1) + rmsFNP_noveto[i - 1];
		}
		else
		{
			y_syst[i - 1] = h_FNP[0]->GetBinContent((i - 1) % NBINSFNP + 1) - rmsFNP_noveto[i - 1];
		}
	}

	g_systematicsRMS_noveto = new TGraph(2 * NBINSFNP, x_syst, y_syst);


	//----------------------------------------------------------------
	// RMS Calculation for FNP With Veto Cut
	//----------------------------------------------------------------

	//Display the systematics as a filled polygon
	float x_syst_veto[NBINSFNP * 2];
	float y_syst_veto[NBINSFNP * 2];

	//Get the RMS of all systematic variations at every pT point
	for (int i = 1; i <= NBINSFNP; i++)
	{
		//Vectors to store, for every FNP bin, the values of the systematic variations that are above and below the default value of FNP
		std::vector<float> valuesAbove;
		std::vector<float> valuesBelow;

		//Array to store, for every FNP bin, the ratio of the systematic variation to the default value
		float valuesRatio[NUMFNPS] = {0.0};
		float values[NUMFNPS] = {0.0};

		for (int j = 0; j < NUMFNPS; j++)
		{
			float fnpDefault = h_FNP_conv[0]->GetBinContent(i);
			float fnp = h_FNP_conv[j]->GetBinContent(i);

			if (fnp > fnpDefault)
			{
				valuesAbove.push_back(fnp - fnpDefault);
			}
			else if (fnp < fnpDefault)
			{
				valuesBelow.push_back(fnp - fnpDefault);
			}

			values[j] = fnp;
			valuesRatio[j] = fnp / h_FNP[0]->GetBinContent(i);
		}

		//Calculate RMS for variations above the default
		rmsFNP_veto_asymm_above[i - 1] = 0.0;
		for (int k = 0; k < valuesAbove.size(); k++)
		{
			rmsFNP_veto_asymm_above[i - 1] += valuesAbove[k] * valuesAbove[k];
		}
		rmsFNP_veto_asymm_above[i - 1] = TMath::Sqrt((float) rmsFNP_noveto_asymm_above[i - 1] / valuesAbove.size());

		//Calculate RMS for variations below the default
		rmsFNP_veto_asymm_below[i - 1] = 0.0;
		for (int k = 0; k < valuesBelow.size(); k++)
		{
			rmsFNP_veto_asymm_below[i - 1] += valuesBelow[k] * valuesBelow[k];
		}
		rmsFNP_veto_asymm_below[i - 1] = TMath::Sqrt((float) rmsFNP_veto_asymm_below[i - 1] / valuesBelow.size());

		float mean = TMath::Mean(NUMFNPS, values);
		float meanRatio = TMath::Mean(NUMFNPS, valuesRatio);
		float rmsFNP_vetoRatio = TMath::RMS(NUMFNPS, valuesRatio);

		b_fnp_syst_veto[i - 1] = new TBox(h_FNP_conv[0]->GetBinCenter(i) - 0.15, h_FNP_conv[0]->GetBinContent(i) - rmsFNP_veto_asymm_below[i - 1], h_FNP_conv[0]->GetBinCenter(i) + 0.15, h_FNP_conv[0]->GetBinContent(i) + rmsFNP_veto_asymm_above[i - 1]);
	}

	//Fill TGraph
	//For a given pT, there will be 2 entries in the TGraph:
	//Central value \pm RMS
	for (int i = 1; i <= 2 * NBINSFNP; i++)
	{
		x_syst_veto[i - 1] = h_FNP[0]->GetBinCenter((i - 1) % NBINSFNP + 1);

		if (i - 1 < NBINSFNP)
		{
			y_syst_veto[i - 1] = h_FNP[0]->GetBinContent((i - 1) % NBINSFNP + 1) + rmsFNP_noveto[i - 1];
		}
		else
		{
			y_syst_veto[i - 1] = h_FNP[0]->GetBinContent((i - 1) % NBINSFNP + 1) - rmsFNP_noveto[i - 1];
		}
	}

	g_systematicsRMS_noveto = new TGraph(2 * NBINSFNP, x_syst_veto, y_syst_veto);
}


void getRMS()
{
	//Display the systematics as a filled polygon
	float x_syst[NBINSFNP * 2];
	float y_syst[NBINSFNP * 2];

	h_rms_plus = (TH1D*) h_FNP[0]->Clone("h_rms_plus");
	h_rms_minus = (TH1D*) h_FNP[0]->Clone("h_rms_minus");

	//Get the RMS of all systematic variations at every pT point
	for (int i = 1; i <= NBINSFNP; i++)
	{
		float values[NUMFNPS] = {0.0};
		float valuesRatio[NUMFNPS] = {0.0};
		for (int j = 0; j < NUMFNPS; j++)
		{
			float fnp = h_FNP[j]->GetBinContent(i);
			values[j] = fnp;
			valuesRatio[j] = fnp / h_FNP[0]->GetBinContent(i);
		}

		rmsFNP_noveto[i - 1] = TMath::RMS(NUMFNPS, values);
		float mean = TMath::Mean(NUMFNPS, values);
		float meanRatio = TMath::Mean(NUMFNPS, valuesRatio);
		float rmsFNP_novetoRatio = TMath::RMS(NUMFNPS, valuesRatio);
		h_rms_plus->SetBinContent(i, meanRatio + rmsFNP_novetoRatio);
		h_rms_minus->SetBinContent(i, meanRatio - rmsFNP_novetoRatio);

		b_fnp_syst_noveto[i - 1] = new TBox(h_FNP[0]->GetBinCenter(i) - 0.15, mean - TMath::RMS(NUMFNPS, values), h_FNP[0]->GetBinCenter(i) + 0.15, mean + TMath::RMS(NUMFNPS, values));
	}

	//Fill TGraph
	//For a given pT, there will be 2 entries in the TGraph:
	//Central value \pm RMS
	for (int i = 1; i <= 2 * NBINSFNP; i++)
	{
		x_syst[i - 1] = h_FNP[0]->GetBinCenter((i - 1) % NBINSFNP + 1);

		if (i - 1 < NBINSFNP)
		{
			y_syst[i - 1] = h_FNP[0]->GetBinContent((i - 1) % NBINSFNP + 1) + rmsFNP_noveto[i - 1];
		}
		else
		{
			y_syst[i - 1] = h_FNP[0]->GetBinContent((i - 1) % NBINSFNP + 1) - rmsFNP_noveto[i - 1];
		}
	}

	g_systematicsRMS_noveto = new TGraph(2 * NBINSFNP, x_syst, y_syst);
}


void getRatiosToDefault()
{
	//Take the ratio of each systematic variation to the default FNP
	for (int i = 0; i < NUMFNPS; i++)
	{
		//Without conversion veto
		h_FNP_default_ratios[i] = (TH1D*) h_FNP[i]->Clone(Form("h_fnp_default_ratio_%i", i));
		h_FNP_default_ratios[i]->Divide(h_FNP[0]);

		h_FNP_default_ratios[i]->SetMarkerColor(kRed);
		h_FNP_default_ratios[i]->SetMarkerStyle(2);
		h_FNP_default_ratios[i]->SetMarkerSize(1.0);

		//With conversion veto
		h_FNP_conv_default_ratios[i] = (TH1D*) h_FNP_conv[i]->Clone(Form("h_fnp_conv_default_ratio_%i", i));
		h_FNP_conv_default_ratios[i]->Divide(h_FNP_conv[0]);

		h_FNP_conv_default_ratios[i]->SetMarkerColor(kBlue);
		h_FNP_conv_default_ratios[i]->SetMarkerStyle(20);
		h_FNP_conv_default_ratios[i]->SetMarkerSize(0.4);
	}

}


void plotFNP(int index)
{
	TCanvas *cFNP = new TCanvas("cFNP", "cFNP", 600, 600);
	formatHistograms(h_FNP[index], "p_{T} [GeV/c]", "F_{NP}", " ");
	formatHistograms(h_FNP_conv[index], "p_{T} [GeV/c]", "F_{NP}", " ");
	formatGraphs(g_eff_FNP[index], "p_{T} [GeV/c]", "F_{NP}", " ");
	formatGraphs(g_eff_FNP_conv[index], "p_{T} [GeV/c]", "F_{NP}", " ");

	g_eff_FNP[index]->SetMarkerStyle(20);
	g_eff_FNP[index]->SetMarkerColor(kAzure + 2);
	g_eff_FNP[index]->SetLineColor(kAzure + 2);
	g_eff_FNP[index]->SetMarkerSize(0.6);

	g_eff_FNP_conv[index]->SetMarkerStyle(20);
	g_eff_FNP_conv[index]->SetMarkerColor(kRed);
	g_eff_FNP_conv[index]->SetLineColor(kRed);
	g_eff_FNP_conv[index]->SetMarkerSize(0.6);

	h_FNP_conv[index]->SetMarkerColor(kRed);
	h_FNP_conv[index]->SetLineColor(kRed);
	h_FNP[index]->GetYaxis()->SetRangeUser(0.0, 1.2);
	h_FNP[index]->GetXaxis()->SetRangeUser(1.0, 10.0);

	g_eff_FNP[index]->GetYaxis()->SetRangeUser(0.0, 1.2);
	g_eff_FNP[index]->GetXaxis()->SetRangeUser(1.0, 10.0);

	if (effErrors)
	{
		g_eff_FNP[index]->Draw("AP");
		g_eff_FNP_conv[index]->Draw("P,same");
	}
	else
	{
		h_FNP[index]->Draw();
		h_FNP_conv[index]->Draw("same");
	}

	//Plot FNP from template fitting method in the first three bins (July 20 - 2017)
	float pT15[3] = {1.25, 1.75, 2.25};

	float fnp_B0[3] = {0.17, 0.25, 0.342};
	float fnp_err_B0[3] = {0.0};

	float fnp_B1[3] = {0.167, 0.26, 0.355};
	float fnp_err_B1[3] = {0.0};

	float err_x[7] = {0.0};

	TGraphErrors *g_fnp_B0 = new TGraphErrors(3, pT15, fnp_B0, err_x, fnp_err_B0);
	TGraphErrors *g_fnp_B1 = new TGraphErrors(3, pT15, fnp_B1, err_x, fnp_err_B1);

	//Plot FNP from template fitting method (May 9 - 2017)
	/*
	float pT15[7] = {1.25, 1.75, 2.25, 2.75, 3.5, 5.0, 7.0};

	float fnp_B0[7] = {0.269985, 0.357814, 0.428488, 0.454788, 0.506252, 0.602844, 0.798613};
	float fnp_err_B0[7] = {0.00361323, 0.00425627, 0.00682485, 0.0111647, 0.014485, 0.0283789, 0.119273};

	float fnp_B1[7] = {0.272355, 0.364761, 0.435008, 0.488402, 0.584216, 0.599103, 0.772028};
	float fnp_err_B1[7] = {0.00250411, 0.0027377, 0.00429937, 0.00675822, 0.00824138, 0.0167099, 0.0551451};

	float err_x[7] = {0.0};

	TGraphErrors *g_fnp_B0 = new TGraphErrors(7, pT15, fnp_B0, err_x, fnp_err_B0);
	TGraphErrors *g_fnp_B1 = new TGraphErrors(7, pT15, fnp_B1, err_x, fnp_err_B1);
	*/

	g_fnp_B0->SetMarkerStyle(20);
	g_fnp_B0->SetMarkerSize(0.7);
	g_fnp_B0->SetMarkerColor(kViolet);
	g_fnp_B0->SetLineColor(kViolet);

	g_fnp_B1->SetMarkerStyle(20);
	g_fnp_B1->SetMarkerSize(0.7);
	g_fnp_B1->SetMarkerColor(kGreen + 2);
	g_fnp_B1->SetLineColor(kGreen + 2);

	g_fnp_B0->Draw("P,same");
	g_fnp_B1->Draw("P,same");

	TLegend *legFNP = new TLegend(0.3, 0.15, 0.85, 0.40);
	legFNP->SetLineColor(kWhite);
	legFNP->AddEntry(h_FNP_conv[index], "F_{NP} With Veto Cut", "PE");
	legFNP->AddEntry(h_FNP[index], "F_{NP} Without Veto Cut", "PE");
	legFNP->AddEntry(g_fnp_B0, "Template Fitting Method in B0 - No Veto Cut", "P");
	legFNP->AddEntry(g_fnp_B1, "F_{NP} From Template Fitting in B1 - No Veto Cut", "PE");
	legFNP->Draw("same");

	//Systematic error band
	//g_systematicsRMS_noveto->SetLineColor(0); //green
	//g_systematicsRMS_noveto->SetLineWidth(0);
	//g_systematicsRMS_noveto->SetFillColorAlpha(kRed, 0.4); //blue
	//g_systematicsRMS_noveto->Draw("f,same"); //draw with Axis and Fill
	//g_systematicsRMS_noveto->Draw("l,same");  //draw only the contour Line
}


void plotFNP()
{
	TCanvas *cFNP = new TCanvas("cFNP", "cFNP", 600, 600);

	for (int index = 0; index < NUMFNPS; index++)
	{
		formatHistograms(h_FNP[index], "p_{T} [GeV/c]", "F_{NP}", " ");
		formatHistograms(h_FNP_conv[index], "p_{T} [GeV/c]", "F_{NP}", " ");
		formatGraphs(g_eff_FNP[index], "p_{T} [GeV/c]", "F_{NP}", " ");
		formatGraphs(g_eff_FNP_conv[index], "p_{T} [GeV/c]", "F_{NP}", " ");

		g_eff_FNP[index]->SetMarkerStyle(20);
		g_eff_FNP[index]->SetMarkerColor(kAzure + 2);
		g_eff_FNP[index]->SetLineColor(kAzure - 9);
		g_eff_FNP[index]->SetMarkerSize(0.6);

		g_eff_FNP_conv[index]->SetMarkerStyle(20);
		g_eff_FNP_conv[index]->SetMarkerColor(kRed);
		g_eff_FNP_conv[index]->SetLineColor(kRed - 9);
		g_eff_FNP_conv[index]->SetMarkerSize(0.6);

		h_FNP_conv[index]->SetMarkerColor(kRed);
		h_FNP_conv[index]->SetLineColor(kRed);
		h_FNP[index]->GetYaxis()->SetRangeUser(0.0, 1.2);
		h_FNP[index]->GetXaxis()->SetRangeUser(1.0, 10.0);

		g_eff_FNP[index]->GetYaxis()->SetRangeUser(0.0, 1.2);
		g_eff_FNP[index]->GetXaxis()->SetRangeUser(1.0, 10.0);

		if (index == 0)
		{
			if (effErrors)
			{
				g_eff_FNP[index]->SetLineColor(kBlue);
				g_eff_FNP_conv[index]->SetLineColor(kRed);

				g_eff_FNP[index]->SetMarkerColor(kBlue);
				g_eff_FNP_conv[index]->SetMarkerColor(kRed);

				g_eff_FNP[index]->Draw("AP");
				g_eff_FNP_conv[index]->Draw("P,same");
			}
			else
			{
				h_FNP[index]->Draw();
				h_FNP_conv[index]->Draw("same");
			}
		}
		else
		{
			if (effErrors)
			{
				g_eff_FNP[index]->Draw("LX,same");
				g_eff_FNP_conv[index]->Draw("LX,same");
			}
			else
			{
				h_FNP[index]->Draw("same");
				h_FNP_conv[index]->Draw("same");
			}
		}
	}

	g_eff_FNP[0]->Draw("P,same");
	g_eff_FNP_conv[0]->Draw("P,same");

	f_fnp_fit_without_veto->SetLineColor(kBlue);
	f_fnp_fit_without_veto->Draw("same");

	f_fnp_fit_with_veto->SetLineColor(kRed);
	f_fnp_fit_with_veto->Draw("same");

	//Plot FNP from template fitting method (May 9 - 2017)
	/*
	float pT15[7] = {1.25, 1.75, 2.25, 2.75, 3.5, 5.0, 7.0};

	float fnp_B0[7] = {0.269985, 0.357814, 0.428488, 0.454788, 0.506252, 0.602844, 0.798613};
	float fnp_err_B0[7] = {0.00361323, 0.00425627, 0.00682485, 0.0111647, 0.014485, 0.0283789, 0.119273};

	float fnp_B1[7] = {0.272355, 0.364761, 0.435008, 0.488402, 0.584216, 0.599103, 0.772028};
	float fnp_err_B1[7] = {0.00250411, 0.0027377, 0.00429937, 0.00675822, 0.00824138, 0.0167099, 0.0551451};

	float err_x[7] = {0.0};

	TGraphErrors *g_fnp_B0 = new TGraphErrors(7, pT15, fnp_B0, err_x, fnp_err_B0);
	TGraphErrors *g_fnp_B1 = new TGraphErrors(7, pT15, fnp_B1, err_x, fnp_err_B1);
	*/

	//Plot FNP from template fitting method in the first three bins (July 20 - 2017)
	float pT15[3] = {1.25, 1.75, 2.25};

	float fnp_B0[3] = {0.17, 0.25, 0.34};
	float fnp_err_B0[3] = {0.0};

	float fnp_B1[3] = {0.173, 0.263, 0.36};
	float fnp_err_B1[3] = {0.0};

	float err_x[7] = {0.0};

	TGraphErrors *g_fnp_B0 = new TGraphErrors(3, pT15, fnp_B0, err_x, fnp_err_B0);
	TGraphErrors *g_fnp_B1 = new TGraphErrors(3, pT15, fnp_B1, err_x, fnp_err_B1);

	g_fnp_B0->SetMarkerStyle(20);
	g_fnp_B0->SetMarkerSize(0.7);
	g_fnp_B0->SetMarkerColor(kViolet);
	g_fnp_B0->SetLineColor(kViolet);

	g_fnp_B1->SetMarkerStyle(20);
	g_fnp_B1->SetMarkerSize(0.7);
	g_fnp_B1->SetMarkerColor(kGreen + 2);
	g_fnp_B1->SetLineColor(kGreen + 2);

	g_fnp_B0->Draw("P,same");
	g_fnp_B1->Draw("P,same");

	for (int i = 0; i < NBINSFNP; i++)
	{
		b_fnp_syst_noveto[i]->Draw("same");
		b_fnp_syst_noveto[i]->SetLineColor(kBlue);
		b_fnp_syst_noveto[i]->SetFillStyle(0);
	}

	for (int i = 0; i < NBINSFNP; i++)
	{
		b_fnp_syst_veto[i]->Draw("same");
		b_fnp_syst_veto[i]->SetLineColor(kRed);
		b_fnp_syst_veto[i]->SetFillStyle(0);
	}

	TLegend *legFNP = new TLegend(0.3, 0.15, 0.85, 0.40);
	legFNP->SetLineColor(kWhite);
	legFNP->AddEntry(g_eff_FNP_conv[0], "F_{NP} (Algebraic) in B0 - With Veto Cut", "P");
	legFNP->AddEntry(g_eff_FNP[0], "F_{NP} (Algebraic) in B1 - No Veto Cut", "P");
	legFNP->AddEntry(g_fnp_B0, "F_{NP} (Template Fitting) in B0 - No Veto Cut", "P");
	legFNP->AddEntry(g_fnp_B1, "F_{NP} (Template Fitting) in B1 - No Veto Cut", "P");
	legFNP->Draw("same");
}


void plotKilledEfficiency()
{
	TCanvas *cKilledEfficiency = new TCanvas("cKilledEfficiency", "cKilledEfficiency", 600, 600);
	formatHistograms(h_killed_efficiency, "p_{T} [GeV/c]", "Survival Rate due to Random Hits", " ");
	h_killed_efficiency->GetYaxis()->SetRangeUser(0.0, 1.0);
	h_killed_efficiency->Draw();
}

/*
void plotSpeciesSurival()
{
	TCanvas *cSurvivalPizeros = new TCanvas("cSurvivalPizeros", "cSurvivalPizeros", 600, 600);
	formatHistograms(h_survival_pizeros, "p_{T} [GeV/c]", "Survival Rate for Electrons from #pi^{0}", " ");
	h_survival_pizeros->GetYaxis()->SetRangeUser(0.0, 1.0);
	h_survival_pizeros->Draw();

	TCanvas *cSurvivalPhotons = new TCanvas("cSurvivalPhotons", "cSurvivalPhotons", 600, 600);
	formatHistograms(h_survival_photons, "p_{T} [GeV/c]", "Survival Rate for Electrons from #gamma", " ");
	h_survival_photons->GetYaxis()->SetRangeUser(0.0, 1.0);
	h_survival_photons->Draw();

	TCanvas *cSurvivalEtas = new TCanvas("cSurvivalEtas", "cSurvivalEtas", 600, 600);
	formatHistograms(h_survival_etas, "p_{T} [GeV/c]", "Survival Rate for Electrons from #eta", " ");
	h_survival_etas->GetYaxis()->SetRangeUser(0.0, 1.0);
	h_survival_etas->Draw();

	if (savePlots)
	{
		cSurvivalPizeros->SaveAs("Plots/SurvivalPizeros.pdf");
		cSurvivalEtas->SaveAs("Plots/SurvivalEtas.pdf");
		cSurvivalPhotons->SaveAs("Plots/SurvivalPhotons.pdf");
	}
}


void plotConversionEfficiency()
{
	TCanvas *cConversionEfficiency = new TCanvas("cConversionEfficiency", "cConversionEfficiency", 600, 600);
	formatHistograms(h_conversion_efficiency, "p_{T} [GeV/c]", "Survival Rate for Photonic Electrons", " ");
	h_conversion_efficiency->GetYaxis()->SetRangeUser(0.0, 1.0);
	h_conversion_efficiency->Draw();

	if (savePlots)
	{
		cConversionEfficiency->SaveAs("Plots/ConversionEfficiency.pdf");
	}
}


void plotDataElectrons()
{
	TCanvas *cElectrons = new TCanvas("cElectrons", "cElectrons", 600, 600);
	cElectrons->SetLogy();
	formatHistograms(h_elec_pT_inclusive, "p_{T} [GeV/c]", "dN_{e}/dp_{T}", " ");
	h_elec_pT_inclusive->Draw();

	TCanvas *cIsolatedElectrons = new TCanvas("cIsolatedElectrons", "cIsolatedElectrons", 600, 600);
	cIsolatedElectrons->SetLogy();
	formatHistograms(h_elec_pT_isolated, "p_{T} [GeV/c]", "dN_{e}^{iso}/dp_{T}", " ");
	h_elec_pT_isolated->Draw();
}


void saveFiles()
{
	TFile *fout = new TFile("fnp_plots.root", "RECREATE");
	h_FNP->Write();
	h_FNP_conv->Write();
	h_conversion_efficiency->Write();
	h_survival_pizeros->Write();
	h_survival_photons->Write();
	h_survival_etas->Write();
}
*/

void plotSystematicRatio()
{
	TCanvas *cRatio = new TCanvas("cRatio", "Ratio - No Veto Cut", 1000, 600);

	formatHistograms(h_FNP_default_ratios[0], "p_{T} [GeV/c]", "Ratio to Default F_{NP}", "");

	h_rms_minus->SetLineColor(kBlack);
	h_rms_minus->SetLineWidth(5);

	h_rms_plus->SetLineColor(kBlack);
	h_rms_plus->SetLineWidth(5);

	h_FNP_default_ratios[0]->Draw("HIST");
	h_FNP_default_ratios[0]->GetYaxis()->SetRangeUser(0.8, 1.2);
	h_FNP_default_ratios[0]->SetMarkerColor(kRed);
	h_FNP_default_ratios[0]->SetLineColor(kRed - 9);
	h_FNP_default_ratios[0]->SetMarkerStyle(20);
	h_FNP_default_ratios[0]->SetMarkerSize(0.4);

	TLegend *legRatio = new TLegend(0.3, 0.15, 0.85, 0.40);
	legRatio->SetLineColor(kWhite);

	legRatio->Draw("same");

	for (int i = 1; i < NUMFNPS; i++)
	{
		h_FNP_default_ratios[i]->Draw("HIST,same");

		//Set various photon variations in different colors
		int photonIndex = getPhotonIndex(i);
		int pionIndex = getPionIndex(i);
		int etaIndex = getEtaIndex(i);

		if (colorCodeVariation == 0)
		{
			if (photonIndex == 0)
			{
				h_FNP_default_ratios[i]->SetLineColor(kBlue);
			}
			else if (photonIndex == 1)
			{
				h_FNP_default_ratios[i]->SetLineColor(kSpring -1);
			}
			else if (photonIndex == 2)
			{
				h_FNP_default_ratios[i]->SetLineColor(kViolet);
			}
			else
			{
				h_FNP_default_ratios[i]->SetLineColor(kOrange - 3);
			}
		}

		if (colorCodeVariation == 1)
		{
			if (pionIndex == 0)
			{
				h_FNP_default_ratios[i]->SetLineColor(kViolet);
			}
			else if (pionIndex == 1)
			{
				h_FNP_default_ratios[i]->SetLineColor(kOrange - 3);
			}
			else if (pionIndex == 2)
			{
				h_FNP_default_ratios[i]->SetLineColor(kSpring - 1);
			}
		}

		if (colorCodeVariation == 2)
		{
			if (etaIndex == 0)
			{
				h_FNP_default_ratios[i]->SetLineColor(kViolet);
			}
			else if (etaIndex == 1)
			{
				h_FNP_default_ratios[i]->SetLineColor(kOrange - 3);
			}
			else if (etaIndex == 2)
			{
				h_FNP_default_ratios[i]->SetLineColor(kSpring - 1);
			}
			else if (etaIndex == 3)
			{
				h_FNP_default_ratios[i]->SetLineColor(kAzure + 7);
			}
			else
			{
				h_FNP_default_ratios[i]->SetLineColor(kPink + 10);
			}
		}
	}

	h_rms_minus->Draw("hist,same");
	h_rms_plus->Draw("hist,same");
}

void fitFNP()
{
	//Fit FNP to a form A(1-exp(-Bx))
	f_fnp_fit_without_veto = new TF1("f_fnp_fit_without_veto", "[0]*(1.0 - TMath::Exp(-1*[1]*x)) + [2]", 1.25, 9.0);
	f_fnp_fit_without_veto->SetParameter(0, 9.25538e-01);
	f_fnp_fit_without_veto->SetParameter(1, 4.53362e-01);
	f_fnp_fit_without_veto->SetParameter(2, -1.98153e-01);
	g_eff_FNP[0]->Fit(f_fnp_fit_without_veto, "Q0R");

	f_fnp_fit_with_veto = new TF1("f_fnp_fit_with_veto", "[0]*(1.0 - TMath::Exp(-1*[1]*x)) + [2]", 1.25, 9.0);
	f_fnp_fit_with_veto->SetParameter(0, 1.33944e+00);
	f_fnp_fit_with_veto->SetParameter(1, 8.02921e-01);
	f_fnp_fit_with_veto->SetParameter(2, -4.59203e-01);
	g_eff_FNP_conv[0]->Fit(f_fnp_fit_with_veto, "R");
	g_eff_FNP_conv[0]->Draw();
}

void calculateSystematicsFNP()
{
	readFiles();
	//setSimErrorsToZero();
	rebinHistograms();
	defineSpectra();
	calculateSpeciesSurvival();
	calculateConversionEfficiency();
	calculateRandomlyKilledEfficiency();
	getCleanElectronSample();
	getFNP();
	getAsymmRMS();
	getRatiosToDefault();
	fitFNP();

	gStyle->SetOptStat(0);
	gStyle->SetErrorX(0);
	//plotSpeciesSurival();
	//plotConversionEfficiency();
	//plotKilledEfficiency();
	//plotDataElectrons();
	plotFNP();
	plotSystematicRatio();
	//saveFiles();
}