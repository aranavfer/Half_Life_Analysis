#include <iostream>
#include <vector>
#include <sstream>
#include <iterator>
#include <string>
#include <stdio.h>
#include <math.h>
#include <TH1.h>
#include <TMath.h>
#include <TF1.h>
#include <TLegend.h>
#include <TCanvas.h>

using namespace std;

void Draw_Comparison(vector<string> Com_Names,vector<double> Com_T12_Wu,vector<double> Com_T12_Err_Wu,vector<double> Com_T12_Ar,vector<double> Com_T12_Err_Ar){
  
  int gordo = Com_Names.size();
  double xxx[gordo];
  double err_x[gordo], w[gordo], we[gordo], a[gordo], ae[gordo];
  string v_name[gordo];
  
  for (int k=0; k<gordo; k++){
    xxx[k] = k;
    err_x[k] = 0.0;
    w[k] = Com_T12_Wu[k];
    we[k] = Com_T12_Err_Wu[k];
    a[k] = Com_T12_Ar[k];
    ae[k] = Com_T12_Err_Ar[k];
    v_name[k] = Com_Names[k];
    
  }
  
  TCanvas *cc = new TCanvas("cc","Half Life Comparison", 200,10,900,700);
  cc->SetGrid();
  
  TGraphErrors *gr_WU = new TGraphErrors(gordo,xxx,w,err_x,we);
  TGraphErrors *gr_AR = new TGraphErrors(gordo,xxx,a,err_x,ae);

      
    for (int i=0;i<gordo;i++){
        int nb = gr_WU->GetXaxis()->FindBin(xxx[i]);
        gr_WU->GetXaxis()->SetBinLabel(nb,Com_Names[i].c_str());
    }
    
    gr_WU->SetTitle("Half Life Comparison");
    gr_WU->SetMarkerColor(4);
    gr_WU->SetMarkerStyle(21);
    gr_WU->SetLineColor(4);
    gr_WU->GetXaxis()->SetTitle("Mass Number");
    gr_WU->GetXaxis()->SetTitleOffset(1.8);
    gr_WU->GetYaxis()->SetTitle("Half Life[s]");
    gr_WU->GetYaxis()->SetRangeUser(-0.5,2.8);
    
    gr_WU->Draw("AP");
    
    gr_AR->SetMarkerColor(8);
    gr_AR->SetMarkerStyle(21);
    gr_AR->SetLineColor(8);
    
    gr_AR->Draw("Psame");
    
    auto legend = new TLegend(0.7,0.7,0.98,0.9);
    legend->AddEntry(gr_WU, "J. Wu","l");
    legend->AddEntry(gr_AR, "All_merged", "l");
    legend->SetTextSize(0.05);
    legend->Draw();
    
    cc->Update();
    
     string uio = "WU_AN_comparison.png";

     
     cc->SaveAs(uio.c_str());
     
     cc->Close();
  
  
}

void Draw_Fit(TH1D* histo1, double X0, double lnda1, double error_l1, double lnda2,double lnda3,double lnda5,double lnda6,double Pb1, double Pb2,double Pb5,double Pn1,double Pn2, double backgrownd, double bta, double bta_err, double chi2){
    
    int Nbin = histo1->GetSize();
    double x_v[Nbin], y1_v[Nbin], y2_v[Nbin], y3_v[Nbin], y5_v[Nbin], y6_v[Nbin], F_v[Nbin], ybac_v[Nbin];
    double bw = histo1->GetBinWidth(0);
    TAxis *xaxis = histo1->GetXaxis();

    TH1D *hres = new TH1D("h1", "", Nbin, histo1->GetXaxis()->GetBinLowEdge(0), histo1->GetXaxis()->GetBinUpEdge(Nbin));
    
    for (int k=0; k<Nbin; k++){
        x_v[k] = 0.0;
        y1_v[k]= 0.0;
        y2_v[k]= 0.0;
        y3_v[k]= 0.0;
        y5_v[k]= 0.0;
        y6_v[k]= 0.0;
        F_v[k] = 0.0;
        ybac_v[k] = backgrownd;
        
        
    }
    
    for (int j=0; j<Nbin; j++){
        double x = xaxis->GetBinCenter(j);
        x_v[j] = x-bw;
        if (x>bw){

            double X1 = X0*exp(-lnda1*(x-bw));
            double X2 = X0*Pb1*lnda1*(exp(-lnda1*(x-bw))/(lnda2-lnda1)+exp(-lnda2*(x-bw))/(lnda1-lnda2));
            double X3 = X0*Pb1*Pb2*lnda1*lnda2*(exp(-lnda1*(x-bw))/((lnda2-lnda1)*(lnda3-lnda1))+exp(-lnda2*(x-bw))/((lnda1-lnda2)*(lnda3-lnda2))+exp(-lnda3*(x-bw))/((lnda1-lnda3)*(lnda2-lnda3)));
    
            double X5 = X0*Pn1*lnda1*(exp(-lnda1*(x-bw))/(lnda5-lnda1)+exp(-lnda5*(x-bw))/(lnda1-lnda5));
            double X6 = X0*Pb1*Pn2*lnda1*lnda2*(exp(-lnda1*(x-bw))/((lnda2-lnda1)*(lnda6-lnda1))+exp(-lnda2*(x-bw))/((lnda1-lnda2)*(lnda6-lnda2))+exp(-lnda6*(x-bw))/((lnda1-lnda3)*(lnda2-lnda3)))+X0*Pn1*Pb5*(exp(-lnda1*(x-bw))/((lnda5-lnda1)*(lnda6-lnda1))+exp(-lnda5*(x-bw))/((lnda1-lnda5)*(lnda6-lnda5))+exp(-lnda6*(x-bw))/((lnda1-lnda6)*(lnda5-lnda6)));

    
             double F = lnda1*X1+lnda2*X2+lnda3*X3+lnda5*X5+lnda6*X6+backgrownd;
            
             y1_v[j] = backgrownd + lnda1*X1;
             y2_v[j] = backgrownd + lnda2*X2;
             y3_v[j] = backgrownd + lnda3*X3;
             y5_v[j] = backgrownd + lnda5*X5;
             y6_v[j] = backgrownd + lnda6*X6;
             F_v[j] = F; 
             
             ybac_v[j] = backgrownd;
	     
 	     hres->SetBinContent(j,histo1->GetBinContent(j)-F_v[j]);
            
        }
        else{
	  
 	  hres->SetBinContent(j,histo1->GetBinContent(j)-backgrownd);
	}
	
    }
    
    TCanvas *c2 = new TCanvas("c2","Half Life Plot",200,10,1000,800);
    TPad *pad1 = new TPad("pad1", "",0.0,0.4,1.0,1.0,21);

    pad1->SetFillColor(0);

    pad1->SetGrid();
    pad1->SetLogy();
    pad1->Draw();
    pad1->cd();
    histo1->Draw();
    histo1->SetTitle(histo1->GetName());
    histo1->GetXaxis()->SetTitle("Time");
    histo1->GetYaxis()->SetTitle("Counts");
    histo1->SetStats(0);
   
    
    TGraph *gr1 = new TGraph(Nbin,x_v,y1_v);
    gr1->SetLineWidth(2);
    gr1->SetLineColor(1);
    gr1->Draw("C");
      
      
    TGraph *gr2 = new TGraph(Nbin,x_v,y2_v);
    gr2->SetLineWidth(2);
    gr2->SetLineColor(8);
    gr2->Draw("C");
      
    TGraph *gr3 = new TGraph(Nbin,x_v,y3_v);
    gr3->SetLineWidth(2);
    gr3->SetLineColor(28);
    gr3->Draw("C");
      
    TGraph *gr5 = new TGraph(Nbin,x_v,y5_v);
    gr5->SetLineWidth(2);
    gr5->SetLineColor(6);
    gr5->Draw("C");
     
    TGraph *gr6 = new TGraph(Nbin,x_v,y6_v);
    gr6->SetLineWidth(2);
    gr6->SetLineColor(46);
    gr6->Draw("C");
    
    TGraph *bac = new TGraph(Nbin,x_v,ybac_v);
    bac->SetLineWidth(2);
    bac->SetLineColor(7);
    bac->Draw("C");
    
       
    TGraph *grf = new TGraph(Nbin,x_v,F_v);
    grf->SetLineWidth(2);
    grf->SetLineColor(2);
    grf->Draw("C");
    
    
    
     string asd = "T12: "+to_string(0.69314718/lnda1)+"+-"+to_string(error_l1)+" [s]";
     string qwe = "Chi2: "+to_string(chi2);
     string zxc = "Backgrownd: "+to_string(backgrownd);
     string jkl = "Beta Eff: "+to_string(bta)+"+-"+to_string(bta_err);
     TPaveText *pt = new TPaveText(4.5,backgrownd+50,10.0,backgrownd+10,"TR");
     pt->AddText(asd.c_str());
     pt->AddText(qwe.c_str());
     pt->AddText(zxc.c_str());
     pt->AddText(jkl.c_str());
     pt->Draw();
     
     
     c2->cd();
     TPad *pad2 = new TPad("pad2", "",0.0,0.0,1.0,0.4,22);
     pad2->SetFillColor(0);
     pad2->SetGrid();
     pad2->Draw();
     pad2->cd();
      
    
//     hres->SetTitle("Residuals");
    hres->GetXaxis()->SetTitle("Time [s]");
    hres->GetYaxis()->SetTitle("Residual");
    hres->GetXaxis()->SetRangeUser(-10,10);
    hres->GetYaxis()->SetRangeUser(-20,20);
    hres->SetStats(0);
    hres->SetLineColor(1);
    hres->Draw("E");
     
     string uio = histo1->GetName();
     uio = uio + ".png";
     
     c2->SaveAs(uio.c_str());
     
     c2->Close();
    
    
    
}

double BateFun(double *x, double *par){
    
    double X1, X2, X3, X5, X6, F;
    double X0, lnda1, lnda2, lnda3, lnda5, lnda6, backgrownd, bw;
    double Pb1, Pb2, Pb5, Pn1, Pn2;
    
    // -------------------------------------------------
    //              Eficiencias
    // -------------------------------------------------
    
    double Beta_Efficiency;
    double ln2 = 0.69314718;
    
    // -------------------------------------------------
    
    X0 = par[0];
    lnda1 = ln2/par[1];
    lnda2 = ln2/par[2];
    lnda3 = ln2/par[3];
    lnda5 = ln2/par[4];
    lnda6 = ln2/par[5];
    
    Pb1 = par[6];
    Pb2 = par[7];
    Pb5 = par[8];
    Pn1 = par[9];
    Pn2 = par[10];
    
    backgrownd = par[11];
    
    bw = par[12];
    
    
    
    X1 = X0*exp(-lnda1*(x[0]-bw));
    X2 = X0*Pb1*lnda1*(exp(-lnda1*(x[0]-bw))/(lnda2-lnda1)+exp(-lnda2*(x[0]-bw))/(lnda1-lnda2));
    X3 = X0*Pb1*Pb2*lnda1*lnda2*(exp(-lnda1*(x[0]-bw))/((lnda2-lnda1)*(lnda3-lnda1))+exp(-lnda2*(x[0]-bw))/((lnda1-lnda2)*(lnda3-lnda2))+exp(-lnda3*(x[0]-bw))/((lnda1-lnda3)*(lnda2-lnda3)));
    
    X5 = X0*Pn1*lnda1*(exp(-lnda1*(x[0]-bw))/(lnda5-lnda1)+exp(-lnda5*(x[0]-bw))/(lnda1-lnda5));
    X6 = X0*Pb1*Pn2*lnda1*lnda2*(exp(-lnda1*(x[0]-bw))/((lnda2-lnda1)*(lnda6-lnda1))+exp(-lnda2*(x[0]-bw))/((lnda1-lnda2)*(lnda6-lnda2))+exp(-lnda6*(x[0]-bw))/((lnda1-lnda3)*(lnda2-lnda3)))+X0*Pn1*Pb5*(exp(-lnda1*(x[0]-bw))/((lnda5-lnda1)*(lnda6-lnda1))+exp(-lnda5*(x[0]-bw))/((lnda1-lnda5)*(lnda6-lnda5))+exp(-lnda6*(x[0]-bw))/((lnda1-lnda6)*(lnda5-lnda6)));

    
    F = lnda1*X1+lnda2*X2+lnda3*X3+lnda5*X5+lnda6*X6+backgrownd;
    
    
    return F;
    
}

vector<double> Half_Life_Fit(TH1D* histo1,int op_back, vector< vector<double> > landa_vec, vector<double> Pb_vec, vector<double> Pn_vec, double Total_implants){

    int Nbin, whichj;
    double fistbin, lastbin, binCenter;
    double Backgrownd, Error_Backgrownd;
    
    TAxis *xaxis = histo1->GetXaxis();
    Nbin = histo1->GetSize();
    
    
    // ------- Buscamos donde se empieza a integrar ---------------


    
    fistbin = xaxis->GetBinCenter(0);
    lastbin = xaxis->GetBinCenter(Nbin);
    
    double bw = histo1->GetBinWidth(0);
    
    for (int j=0; j<Nbin; j++){
        binCenter = xaxis->GetBinCenter(j);
        if (0 < binCenter && binCenter < bw){
            whichj = j;

            
        }
        
    }
   
    
    // --------------- CALCULO BACKGROWND ------------------------
    
    if (op_back == 0){
        histo1->Fit("pol0","Lq0","0",fistbin,0);
        Backgrownd = histo1->GetFunction("pol0")->GetParameter(0);
        Error_Backgrownd = histo1->GetFunction("pol0")->GetParError(0);
    }
    else if (op_back == 1){
        cout << "Opcion Backgrownd No disponible todavia " << endl;
    }
    

    // ----------------- HISTOGRAM FIT ----------------------------
    
    TF1 *f1 = new TF1("myfunc", BateFun, 0.0, lastbin, 13);
    
    double X0,X00, lnda1, lnda2, lnda3, lnda5, lnda6, Pb1, Pb2, Pb5, Pn1, Pn2;
    double err_low_lnda1, err_low_lnda2, err_low_lnda3, err_low_lnda5, err_low_lnda6;
    double err_up_lnda1, err_up_lnda2, err_up_lnda3, err_up_lnda5, err_up_lnda6;
    double Pb1_val, Pb2_val, Pb5_val, Pn1_val, Pn2_val;
    
    X00 = histo1->GetBinContent(whichj+1);
    
    // -------------- Las variables se llaman landa, pero en realidad son T12 ------------------ 
    
    lnda1 = landa_vec[0][0]; err_low_lnda1 = landa_vec[0][1]; err_up_lnda1 = landa_vec[0][2];
    lnda2 = landa_vec[1][0]; err_low_lnda2 = landa_vec[1][1]; err_up_lnda2 = landa_vec[1][2];
    lnda3 = landa_vec[2][0]; err_low_lnda3 = landa_vec[2][1]; err_up_lnda3 = landa_vec[2][2];
    lnda5 = landa_vec[3][0]; err_low_lnda5 = landa_vec[3][1]; err_up_lnda5 = landa_vec[3][2];
    lnda6 = landa_vec[4][0]; err_low_lnda6 = landa_vec[4][1]; err_up_lnda6 = landa_vec[4][2];
    
    Pb1_val = Pb_vec[0]; Pb2_val = Pb_vec[1]; Pb5_val = Pb_vec[3];
    Pn1_val = Pn_vec[0]; Pn2_val = Pn_vec[1];

    f1->SetParameter(0,X00);
    f1->SetParameter(1, lnda1);
    f1->FixParameter(6, Pb1_val);
    f1->FixParameter(7, Pb2_val);
    f1->FixParameter(8, Pb5_val);
    f1->FixParameter(9, Pn1_val);
    f1->FixParameter(10, Pn2_val);
    f1->FixParameter(12, bw);
    
    // --------------------------------------------------
    
    f1->SetParName(0,"X0");
    f1->SetParName(1,"T12_P");
    f1->SetParName(2,"T12_D");
    f1->SetParName(3,"T12_G_D");
    f1->SetParName(4, "T12_n_D");
    f1->SetParName(5, "T12_n_G_D");
    f1->SetParName(6, "Pb1");
    f1->SetParName(7, "Pb2");
    f1->SetParName(8, "Pb5");
    f1->SetParName(9, "Pn1");
    f1->SetParName(10, "Pn2");
    f1->SetParName(11, "Backgro");
    f1->SetParName(12,"Bin_width");

    
    // lnda2, lnda3, lnda5, lnda6 are going to be gaussian variables.
    TRandom3 *r2 = new TRandom3();
    r2->SetSeed(0);
    TRandom3 *r3 = new TRandom3();
    r3->SetSeed(0);
    TRandom3 *r5 = new TRandom3();
    r5->SetSeed(0);
    TRandom3 *r6 = new TRandom3();
    r6->SetSeed(0);
    TRandom3 *r7 = new TRandom3();
    r7->SetSeed(0);

    
    int Nrolls = 1000;
    vector<double> v_t12_fit, v_error_fit, v_chichi2;
    
    double mean_t12 = 0; 
    double mean_t12_error = 0;
    double sigma_fit = 0;
    double mean_chichi = 0;
    double mean_X0_error = 0;
    

    
    for (int rr=0; rr<Nrolls; rr++){
        double ra_2 = (r2->Gaus(lnda2, err_low_lnda2));
        double ra_3 = (r3->Gaus(lnda3, err_low_lnda3));
        double ra_5 = (r5->Gaus(lnda5, err_low_lnda5));
        double ra_6 = (r6->Gaus(lnda6, err_low_lnda6));
        double ra_7 = (r7->Gaus(Backgrownd, Error_Backgrownd));
        
        
        if (ra_2 <= 0.0) continue;
        if (ra_3 <= 0.0) continue;
        if (ra_5 <= 0.0) continue;
        if (ra_6 <= 0.0) continue;
        
        f1->FixParameter(2, ra_2);
        f1->FixParameter(3, ra_3);
        f1->FixParameter(4, ra_5);
        f1->FixParameter(5, ra_6);
        f1->FixParameter(11, ra_7);

        
        histo1->Fit("myfunc","LR+0","",bw,lastbin);
        
	double error_fit, chichi2;
	
        double t12_fit = f1->GetParameter(1);
        double X0_fit = f1->GetParameter(0);
	double X0_error = f1->GetParError(0);
        

        error_fit = f1->GetParError(1);
     
        chichi2 = f1->GetChisquare();
        
  
        v_t12_fit.push_back(t12_fit);
        v_error_fit.push_back(error_fit);
        v_chichi2.push_back(chichi2);
        
        mean_t12 = mean_t12+t12_fit/Nrolls;
	mean_t12_error = mean_t12_error+ error_fit/Nrolls;
	mean_X0_error = mean_X0_error + X0_error/Nrolls;
        mean_chichi = mean_chichi+chichi2/Nrolls;
        X0 = X0 + X0_fit/Nrolls;

        
    }
    
    for (int rk=0; rk<Nrolls; rk++){
        sigma_fit = sigma_fit + (v_t12_fit[rk]-mean_t12)*(v_t12_fit[rk]-mean_t12)/Nrolls;
    }
    

    sigma_fit = sqrt(sigma_fit);
    
    sigma_fit = sqrt(sigma_fit*sigma_fit+mean_t12_error*mean_t12_error);
    
    vector<double> Function_Output;
    
    // ------------------ CALCULO BETA EFFICIENCY --------------- //
    

    double mean_bta = X0/(Total_implants*xaxis->GetBinWidth(0));

    double sigma_bta = 0;
    
   
    
    Function_Output.push_back(X0);
    Function_Output.push_back(mean_t12);
    Function_Output.push_back(sigma_fit);
    Function_Output.push_back(Backgrownd);
    Function_Output.push_back(mean_chichi);
    Function_Output.push_back(mean_bta);
    Function_Output.push_back(sigma_bta);
    
    return Function_Output;

}

vector<double> Pn_Fit(TH1D* histo1, TH1D* histo2, double bta_eff){

  
  int Nbin = histo2->GetSize();

  double xfini = histo2->GetXaxis()->FindBin(0.0);
  double Back_Pn = histo2->Integral(0,xfini-1);
  double Tnt_Pn = histo2->Integral(xfini+1,Nbin)-Back_Pn;
  
  int Nbin2 = histo1->GetSize();
  double xfini2 = histo1->GetXaxis()->FindBin(0.0);
  double Back_Bt = histo1->Integral(0,xfini2-1);
  double Tnt_Bt = histo1->Integral(xfini2+1,Nbin2)-Back_Bt;
  
  cout << Tnt_Pn << endl;
  cout << Tnt_Bt << endl;
  
  double Pn_value = Tnt_Pn/Tnt_Bt*(0.62/bta_eff);
  
  vector<double> PN_sol;

  PN_sol.push_back(Pn_value);

  return PN_sol;
  
  
  
  
}

void FITTING_PROGRAM(string file_Name){
    
    string dummy;
    
    // ----------- Read Half Life Data Base ----------------
    
    string Data_File_Name = "Half_Life_Data";
    
    ifstream Half_Life_Data(Data_File_Name.c_str());
    
    getline(Half_Life_Data, dummy);
    getline(Half_Life_Data, dummy);
    getline(Half_Life_Data, dummy);
    getline(Half_Life_Data, dummy);
    
    int Number_Nuclei;
    
    Half_Life_Data >> dummy >> Number_Nuclei;
    
    string Data_Names[Number_Nuclei];
    double Data_T12[Number_Nuclei], Data_Error_Low[Number_Nuclei], Data_Error_Upper[Number_Nuclei];
    double Data_Prob_Bta[Number_Nuclei], Data_Prob_Ntr[Number_Nuclei];
    
    string name;
    double t12, err1, err2, pb, pn;
    
    for (int i=0; i<Number_Nuclei; i++){
        Half_Life_Data >> name >> t12 >> err1 >> err2 >> pb >> pn;
        Data_Names[i] = name;
        Data_T12[i] = t12;
        Data_Error_Low[i] = err1;
        Data_Error_Upper[i] = err2;
        Data_Prob_Bta[i] = pb;
        Data_Prob_Ntr[i] = pn;
        
    }
    
    Half_Life_Data.close();
    
    // ------------ Read Reaction Chain Data Base ----------
    
    string Chain_File_Name = "Chain_Reaction_Data";
    ifstream Chain_Life_Data(Chain_File_Name.c_str());
    
    getline(Chain_Life_Data, dummy);
    getline(Chain_Life_Data, dummy);
    getline(Chain_Life_Data, dummy);
    getline(Chain_Life_Data, dummy);
    
    int Number_Reaction, Number_Decays;
    
    Chain_Life_Data >> dummy >> Number_Reaction >> dummy >> Number_Decays;
    
    string Chain_Reactions[Number_Reaction][Number_Decays+1];
    
    for (int j=0; j<Number_Reaction; j++){
        for (int k=0; k<Number_Decays+1; k++){
            Chain_Life_Data >> name;
            Chain_Reactions[j][k] = name;
        }
    }
    Chain_Life_Data.close();
  
    // ------------- Read Initial Parameters File -----------
    
    
    ifstream Initial_Parameters_File(file_Name.c_str());
    
    string Parent_Name;
    string Calculo_T12, Calculo_Pn;
    int Opcion_Backgrownd, N_of_Analysis;
    string name_data_file, name_data_file1, name_data_file2;
    
    getline(Initial_Parameters_File, dummy);
    getline(Initial_Parameters_File, dummy);
    getline(Initial_Parameters_File, dummy);
    getline(Initial_Parameters_File, dummy);
    getline(Initial_Parameters_File, dummy);
    Initial_Parameters_File >> dummy >> dummy >> dummy >>  N_of_Analysis;
    
    vector <string> V_Names;
    getline(Initial_Parameters_File, dummy);
    for (int k=0; k<N_of_Analysis; k++){
	string name;
	getline(Initial_Parameters_File, name);
	V_Names.push_back(name);
	
    }

    getline(Initial_Parameters_File, dummy);

    getline(Initial_Parameters_File, Calculo_T12);
    getline(Initial_Parameters_File, dummy);
   
    getline(Initial_Parameters_File, Calculo_Pn);
    getline(Initial_Parameters_File, dummy);

    Initial_Parameters_File >> Opcion_Backgrownd;
    getline(Initial_Parameters_File, dummy);
    getline(Initial_Parameters_File,dummy);
    Initial_Parameters_File >> name_data_file;
  
    
    name_data_file1 = name_data_file+".log";
    name_data_file2 = name_data_file+".root";
    
    cout << name_data_file2 << endl;
    TFile *file0 = TFile::Open(name_data_file2.c_str());
    
    // -------------- Total number of implants ------------ //
    ofstream outfile("OUTPUT.txt");
    outfile << "Summary of the results from the T12 and Pn calculations" << endl;
    outfile << "    " << endl;
    outfile << "Nuclei" << "  |  " << "Background" << "  |  " << "Chi2" << "  |  " <<  "Half Life [J. Wu]" << "  |  ";
    outfile << "Half Life [Analysis]" << "  |  " << "Beta Eff." << "  |  " << "Pn [Möller]" << "  |  " << "Pn [Analysis]" << endl;
    outfile << " -------------------------------------------------------------------------------------------------------" ;
    outfile << "--------------------------------------" << endl;

    vector<double> Com_T12_Wu, Com_T12_Err_Wu, Com_T12_Ar, Com_T12_Err_Ar;
    vector<string> Com_Names;
    
    for (int kkkk=0; kkkk<N_of_Analysis; kkkk++){
      
     Parent_Name = V_Names[kkkk];
      
    
    
    ifstream implantsfile(name_data_file1.c_str());
    getline(implantsfile, dummy);
    int N_implant = 0;
    string imname; 
    int imp;
    int pid;
    
  
    for (int i=0; i<38; i++){
      implantsfile >> imname >> imp >> pid;
      if (Parent_Name == imname){
	N_implant = imp;
      }
    }
    
    
    // ------------- Which is our reaction ---------------
    int which_reaction = 999999999;
    
    for (int i=0; i<Number_Reaction; i++){
        if (Parent_Name == Chain_Reactions[i][0]){
            which_reaction = i;
        }
    }
    
    // ----------------------------------------------------
    
    if (which_reaction == 999999999){
        cout << "------------------------------------------------" << endl;
        cout << Parent_Name << " Is not in the Chain_Reaction_Data File" << endl;
        cout << "------------------------------------------------" << endl;
    }
    
    // ----------------------------------------------------------------- \\
    //              Obtencion datos necesarios                           \\
    // ----------------------------------------------------------------- \\
    // La funcion de fit T12 necesita como input:                        \\
    // Nombre del histograma Parent_Name                                 \\
    // Matriz con parametros necesarios [nuclei][t12 err1 err2 Pb Pn]    \\
    // Otros Parametros [Eff_medida_Beta  Eff_medida_Neutrones]          \\
    // ----------------------------------------------------------------- \\
    
    TCanvas *ci = new TCanvas("c2","Implantation Plot",200,10,1000,800);
    string imp_name = "hZi"+Parent_Name;
    TH1D *Himp = (TH1D*)gDirectory->Get(imp_name.c_str());
    Himp->Draw();
    string rty = Himp->GetName();
     rty = rty + ".png";
     
     ci->SaveAs(rty.c_str());
     ci->Close();
    
    
    string Histogram_Name = "hTib"+Parent_Name;
    string Histogram_Pn_Name = "hTibn"+Parent_Name;
    TH1D *Histo = (TH1D*)gDirectory->Get(Histogram_Name.c_str());
    TH1D *Histo_Pn = (TH1D*)gDirectory->Get(Histogram_Pn_Name.c_str());
    
    
     Histo->Rebin(6);
    
    vector< vector<double> > Vector_of_Lndas_and_errors;
    vector<double> Vector_of_Pbs;
    vector<double> Vector_of_Pns;
    vector<double> v_each_nuclei;
    
    cout << endl;
    cout << "--------------------------------------------------------" << endl;
    cout << "The reaction Chain of " << Parent_Name << " is: " << endl;
    for (int k=0; k<Number_Decays+1; k++){
        cout << Chain_Reactions[which_reaction][k] << "      ";
    }
    cout << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << endl;
    
    
    for (int k=0; k<Number_Decays+1; k++){
        int which_pos = 999999999;
        
        for (int j=0; j<Number_Nuclei; j++){
            
            if (Chain_Reactions[which_reaction][k] == Data_Names[j]){
                which_pos = j;
            }
        }
        
        if  (which_pos == 999999999){
            cout << "--------------------------------------" << endl;
            cout << Chain_Reactions[which_reaction][k] << " Is not in the Half_Life_Data File" << endl;
            cout << "---------------------------------------" << endl;
            break;
        }
        
        v_each_nuclei.push_back(Data_T12[which_pos]);
        v_each_nuclei.push_back(Data_Error_Low[which_pos]);
        v_each_nuclei.push_back(Data_Error_Upper[which_pos]);
        
        Vector_of_Pbs.push_back(Data_Prob_Bta[which_pos]);
        Vector_of_Pns.push_back(Data_Prob_Ntr[which_pos]);
        
        Vector_of_Lndas_and_errors.push_back(v_each_nuclei);
        v_each_nuclei.clear();
        
    }
    
    // ---------------------------- CALCULATION T12 and PN --- OUTPUT GENERATION ---------------------------------- //
    
    
    double bta_out = 0;
    
    if (Calculo_T12 == "1"){

           vector<double> T12 = Half_Life_Fit(Histo, Opcion_Backgrownd, Vector_of_Lndas_and_errors, Vector_of_Pbs, Vector_of_Pns, N_implant);
           
           double t12_out, error_out, Back_out, chi_out, bta_error, x0;
           
           x0 = T12[0];
           t12_out = T12[1]; error_out = T12[2];
           Back_out = T12[3]; chi_out = T12[4];
           bta_out = T12[5]; bta_error = T12[6];
           
            cout << " ------------------ Half Life --------------------- " << endl;
            cout << t12_out << "+-" << error_out << endl;
            cout << "---------------------------------------------------" << endl;
            cout << "Backgrownd: " << Back_out << "      Chi2: " << chi_out << endl;
            cout << "----------------------------------------------------" << endl;
            cout << "Beta Efficiency: "<< bta_out << "+-" << bta_error << endl;
            cout << "----------------------------------------------------" << endl;
            
            
            double lnda1_plot = 0.69314718/t12_out;
            double lnda2_plot = 0.69314718/Vector_of_Lndas_and_errors[1][0]; 
            double lnda3_plot = 0.69314718/Vector_of_Lndas_and_errors[2][0];
            double lnda5_plot = 0.69314718/Vector_of_Lndas_and_errors[3][0];
            double lnda6_plot = 0.69314718/Vector_of_Lndas_and_errors[4][0];
            double Pb1_plot = Vector_of_Pbs[0]; 
            double Pb2_plot = Vector_of_Pbs[1]; 
            double Pb5_plot = Vector_of_Pbs[3];
            double Pn1_plot = Vector_of_Pns[0]; 
            double Pn2_plot = Vector_of_Pns[1];
            
     
            
             Draw_Fit(Histo, x0, lnda1_plot, error_out, lnda2_plot, lnda3_plot, lnda5_plot, lnda6_plot, Pb1_plot, Pb2_plot, Pb5_plot, Pn1_plot, Pn2_plot, Back_out, bta_out, bta_error, chi_out);
            
	     double t12_wu = Vector_of_Lndas_and_errors[0][0];
	     double t12_err_wu = max(Vector_of_Lndas_and_errors[0][1],Vector_of_Lndas_and_errors[0][2]);
	     double pn_moller = Vector_of_Pns[0];
	     
	    outfile << Parent_Name << "      " << Back_out << "       " << chi_out << "       ";
	    outfile << t12_wu << "+-" << t12_err_wu << "      " << t12_out << "+-" << error_out;
	    outfile << "      " << bta_out << "+-" << bta_error << "       ";
	    outfile << pn_moller << "       ";
	    
	    Com_T12_Wu.push_back(t12_wu);
	    Com_T12_Err_Wu.push_back(t12_err_wu);
	    Com_T12_Ar.push_back(t12_out);
	    Com_T12_Err_Ar.push_back(error_out);
	    Com_Names.push_back(Parent_Name);
	    
	    
    }

    if (Calculo_Pn == "1"){
      vector<double> PPPn = Pn_Fit(Histo,Histo_Pn, bta_out);

      cout << "-------------------- PN VALUE -----------------" << endl;
      cout << PPPn[0] << endl;
      cout << "-----------------------------------------------" << endl;
      
      outfile << PPPn[0] << endl;
      
    }
    else{
      outfile << "-" << endl;
    }
    
    
    
    
    
    // ---------------- PREPARATION FOR NEXT NUCLEI -------------------------- //
    
    Vector_of_Pns.clear();
    Vector_of_Pbs.clear();
    Vector_of_Lndas_and_errors.clear();
    v_each_nuclei.clear();
    
    }
    outfile.close();
    Draw_Comparison(Com_Names,Com_T12_Wu, Com_T12_Err_Wu, Com_T12_Ar, Com_T12_Err_Ar);
    
}
