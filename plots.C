#ifdef __MAKECINT__
#pragma link C++ class vector<vector<int> >+;
#pragma link C++ class vector<vector<float> >+;
#endif

#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <iomanip>
#include <cmath>
#include <fstream>

// Fonction pour trier un vecteur de paires en fonction du premier paramètre
bool myComparison(const pair<double,double> &a,const pair<double,double> &b)
{
    return a.first<b.first;
}

//
// Fonction utilisant les données du vecteur data pour calculer l'evolution du paramètre au cours du temps
// La période d'étude est définie par tmin et tmax (temps en secondes)
// Les données sont moyennées sur une période de ndays jours glissants, correspondant à delay secondes
// Les résultats normalisés ou pas sont stockés dans les histogrammes ahisto et ahisto_raw respectivement
//

void doPlot(const std::vector<pair<double,double> > *data,TH2F *ahisto,TH2F *ahisto_raw, double tmin, double tmax, int ndays)
{
    double sum   = 0;
    double mtime = 0;
    double nent = 0;
    double date,SC,mult;
    double delay=ndays*24*3600;
    
    int npts= ndays*(int((tmax-tmin)/delay)+1); // Nombre de points de mesure sur la période (1/jour)
    
    std::vector<double> time;      // Nombre de points sur un jour (peu être >1 car il y a plusieurs schémas vaccinaux)
    std::vector<double> values;    // Valeur brute sur la période
    std::vector<double> pop;       // Population totale représentée
    time.clear();
    values.clear();
    pop.clear();
    
    for(int i = 0; i < npts; ++i)
    {
        time.push_back(0);
        values.push_back(0);
        pop.push_back(0);
    }
        
    int idx;
    
    // Boucle sur les données, on remplit les
    // trois vecteurs
    
    for(int i = 0; i < int(data->size()); ++i)
    {
        date=data->at(i).first;
        idx=int((date-tmin)/(delay/ndays));
        mult=long(data->at(i).second/1000.);
        SC=data->at(i).second-mult*1000;
        if (SC>=999) SC=SC-1000.; // Minor rouding bug I suppose
            
        mult=float(mult)/1e9;

        time[idx]+=1;
        values[idx]+=SC*mult; // On dénormalise la valeur mesurée
        pop[idx]+=mult;
    }

    // On a terminé, on calcule les moyennes glissantes
    for(int i = 0; i < npts-ndays; ++i)
    {
        mult=0;
        sum=0;
        nent=0;
        mtime=(tmin+(float(i)+0.5)*(delay/ndays));
        
        // Calcul de la moyenne glissante sur ndays jours
        for(int j = 0; j < ndays; ++j)
        {
            // On vérifie qu'il y a des données pour ce jour...
            if (time[i+j]==0) continue;
            
            mult+=values[i+j];         // Somme totale
            sum+=values[i+j]/pop[i+j]; // Somme normalisée (par million)
            nent+=1;
        }
        
        
        if (mtime<ndays-2) continue; // Si il manque plus de 2 jours sur la fenêtre on passe
        
        ahisto->Fill(mtime,sum/nent);      // Moyenne totale
        ahisto_raw->Fill(mtime,mult/nent); // Moyenne des valeurs normalisées
    }
}



//
// Fonction utilisant les données des vecteur data (vax) et datan (non-vax) pour calculer l'evolution
// du rapport non-vax/vax au cours du temps, pour un age donné
// La période d'étude est définie par tmin et tmax (temps en secondes)
// Les données sont moyennées sur une période de ndays jours glissants, correspondant à delay secondes
// Les résultats sont stockés dans les histogrammes histo et histo2
//

void doRatioPlot(const std::vector<pair<double,double> > *data,const std::vector<pair<double,double> > *datan,
                 TH2F *histo,TH2F *historaw,int age,TH2F *histo2,TH2F *histo2raw, double tmin, double tmax,int ndays)
{
    double sum   = 0;
    double mtime = 0;
    double sum_nvac   = 0;
    double mult_nvac = 0;
    
    double date,SC,mult;
    double delay=ndays*24*3600;
    
    int npts= ndays*(int((tmax-tmin)/delay)+1); // Nombre de points de mesure sur la période (1/jour)

    
    std::vector<double> time;         // Nombre de points sur un jour (peu être >1 car il y a plusieurs schémas vaccinaux)
    std::vector<double> values;       // Valeur brute sur la période
    std::vector<double> pop;          // Population totale représentée
    std::vector<double> time_nvac;    //
    std::vector<double> values_nvac;  // Idem pour les non-vaccinées
    std::vector<double> pop_nvac;     //
    
    time.clear();
    values.clear();
    pop.clear();
    time_nvac.clear();
    values_nvac.clear();
    pop_nvac.clear();

    for(int  i = 0; i < npts; ++i)
    {
        time.push_back(0);
        values.push_back(0);
        pop.push_back(0);
        time_nvac.push_back(0);
        values_nvac.push_back(0);
        pop_nvac.push_back(0);
    }
    
    int idx;
    
    // On commence par collecter les données pour les vaccinés
    for(int  i = 0; i < int(data->size()); ++i)
    {
        date=data->at(i).first;
        idx=int((date-tmin)/(delay/ndays));
        mult=long(data->at(i).second/1000.);
        SC=data->at(i).second-mult*1000;
        if (SC>=999) SC=SC-1000.; // Minor rounding bug I suppose
            
        mult=float(mult)/1e9;

        time[idx]+=1;
        values[idx]+=SC*mult;
        pop[idx]+=mult;
    }
        
    // Ensuite les données pour les non-vaccinés
    for(int i = 0; i < int(datan->size()); ++i)
    {
        date=datan->at(i).first;
        idx=int((date-tmin)/(delay/ndays));
        mult=long(datan->at(i).second/1000.);
        SC=datan->at(i).second-mult*1000;
        if (SC>=999) SC=SC-1000.; // Minor rouding bug I suppose
                
        mult=float(mult)/1e9;
        
        time_nvac[idx]+=1;
        values_nvac[idx]+=SC*mult; // Weighted average (SC is already in per millions
        pop_nvac[idx]+=mult;
    }
            
    double vac,nonvac,vacraw,nonvacraw,rapport;
    
    // Puis on calcule les moyennes glissantes
    for(int i = 0; i < npts-ndays; ++i)
    {
        rapport=0;
        vac=0;
        nonvac=0;
        vacraw=0;
        nonvacraw=0;
        mult=0;
        mtime=(tmin+(float(i)+0.5)*(delay/ndays));
        
        // Moyenne glissante sur ndays
        for(int j = 0; j < ndays; ++j)
        {
            // Evidemment si il manque des données pour un point on arrête
            if (time_nvac[i+j]==0 || time[i+j]==0) continue;
            
            nonvac+=values_nvac[i+j]/pop_nvac[i+j]; // Non-vaccinés, entrées par millions
            vac+=values[i+j]/pop[i+j];              // Vaccinés, entrées par millions
            nonvacraw+=values_nvac[i+j];            // Non-vaccinés
            vacraw+=values[i+j];                    // Vaccinés
            mult+=1;
        }
        if (mult<ndays-2) continue; // Si il manque plus de 2 jours sur la fenêtre on passe
            
        // Remplissage des histos
        histo->Fill(mtime,nonvac/vac);
        historaw->Fill(mtime,nonvacraw/vacraw);
        histo2->Fill(mtime,age,nonvac/vac);
        histo2raw->Fill(mtime,age,nonvacraw/vacraw);
    }
}

//
// Fonction générant un fichier root contenant toutes les données utiles
// à partir du fichier texte produit par la macro python de parsing
//

void write(const TString& filename)
{
    TFile hfile("vaxcovid.root","RECREATE","");
    TTree *tree = new TTree("DreesData","");

    double date=0;
    double age=0;
    double vac=0;
    double HO=0;
    double HOp=0;
    double SC=0;
    double SCp=0;
    double pop=0;
    
    tree->Branch("date_o",&date,"date_o/D");  // Date (en secondes)
    tree->Branch("age_o",&age,"age_o/D");     // Classe d'age
    tree->Branch("vac_o",&vac,"vac_o/D");     // Statut vaccinal
    tree->Branch("HO_o",&HO,"HO_o/D");        // Hospitalisation
    tree->Branch("HOp_o",&HOp,"HOp_o/D");     // Hospitalisation avec PCR+
    tree->Branch("SC_o",&SC,"SC_o/D");        // Soins critiques
    tree->Branch("SCp_o",&SCp,"SCp_o/D");     // Soins critiques avec PCR+
    tree->Branch("pop_o",&pop,"pop_o/D");     // Population dans la catégorie (en millions)

    /*
    Categories (45 catégories)
     
    --> Statut vaccinal (vac)
     
     0 = Complet de 6 mois et plus - avec rappel
     1 = Complet de 6 mois et plus - sans rappel
     2 = Complet de moins de 3 mois - avec rappel
     3 = Complet de moins de 3 mois - sans rappel
     4 = Complet entre 3 mois et 6 mois - avec rappel
     5 = Complet entre 3 mois et 6 mois - sans rappel
     6 = Non-vaccinés
     7 = Primo dose efficace  // Pas consideré
     8 = Primo dose récente   // Pas consideré
     
     --> Tranche d'age
     
     0 = [40,59]
     1 = [80+]
     2 = [0,19]               // Pas consideré
     3 = [20,39]
     4 = [60,79]
     */
    
        
    ifstream in(filename);
    if (!in) return;
    
    int ages=5;
    int evtlength=10;
    int compt=0;
    std::string line;
    
    // Read the text file to retrieve the parameters
    while (std::getline(in, line))
    {
        if (compt==0)
        {
            in >> date;
            ++compt;
        }
        else
        {
            age=int((compt-1)/evtlength);
            vac=int((compt-evtlength*age-1));
            in >> HO;
            in >> HOp;
            in >> SC;
            in >> SCp;
            in >> pop;
            ++compt;
            
            tree->Fill();
        }
        if (compt==(evtlength*ages)+1) compt=0;
    }
    
    tree->Write();
    in.close();
    hfile.Write();
    hfile.Close();
     
}


//
// Fonction réalisant l'analyse des données
//
// filename est le fichier root produit par la fonction write
// minbycat est le nombre minimum de personnes dans la catégorie pour prendre en compte le point (en million)
// ndays est la taille de la fenêtre sur laquelle la moyenne est calculée (en jours). O
//
// 4 variables sont disponible:
// -> HO  (0) : nombre de personnes hospitalisées
// -> HOp (1) : nombre de personnes hospitalisées avec test PCR+
// -> SC  (2) : nombre de personnes en soins critiques
// -> SCp (3) : nombre de personnes en soins critiques avec test PCR+
//
//
// On calcule ensuite un moyenne glissante pour chaque jour des parametres suivants:
//
// Variable normalisée en fonction du temps (par millions d'habitants): vaccinés et non-vaccinés
// Variable brute en fonction du temps: vaccinés et non-vaccinés
// Rapport var normalisées en fonction du temps (par millions d'habitants): non-vax/vax
// Rapport var brutes en fonction du temps: non-vax/vax
//
//

void read(std::string filename, float minbycat, int ndays, int var)
{
  if (var<0 || var>3)
  {
      std::cout << "Var should be between 0 and 3 !!!";
      return;
  }
    
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
    
  double date;
  double age;
  double vac;
  double pop;
  double val=0.;

  double delay=ndays*24*3600; // Intervalle de temps en secondes
            
  TChain *newtree = new TChain("DreesData");
  newtree->Add(filename.c_str());
  newtree->SetBranchAddress("date_o",&date);
  newtree->SetBranchAddress("age_o",&age);
  newtree->SetBranchAddress("vac_o",&vac);
  if (var==0) newtree->SetBranchAddress("HO_o",&val);
  if (var==1) newtree->SetBranchAddress("HOp_o",&val);
  if (var==2) newtree->SetBranchAddress("SC_o",&val);
  if (var==3) newtree->SetBranchAddress("SCp_o",&val);
  newtree->SetBranchAddress("pop_o",&pop);

  int nentries=newtree->GetEntries();
  
  double tmin     = 1e20;
  double tmax     = -1;
  double tmaxreal = -1;
  int npts        = 0;
  int npts_r      = 0;
  double valmax   = 0;
  double rawvalmax= 0;
    
  // Seuil minimum pour cette variable
  if (minbycat<=0.01) minbycat=0.01;
    
  // On fait un premier tour pour déterminer les limites
  for (int i=0;i<nentries;++i)
  {
      newtree->GetEntry(i);
      if (vac>6) continue;        // Regarde que complet ou non-vacciné
      if (pop<minbycat) continue; // Pas assez de données (<500k)
      if (age==2) continue;       // Catégorie pas étudiée car peu vaccinée
      
      if (date<tmin)         tmin=date;
      if (date>tmax)         tmax=date;
      if (val>valmax)        valmax=val;
      if (val*pop>rawvalmax) rawvalmax=val*pop;
  }
        
  npts     = int((tmax-tmin)/delay);
  tmaxreal = tmax;
  tmax     = tmin+(npts+1)*delay;
  npts_r   = (tmaxreal-tmin)/(delay/ndays)+1;
    
  std::cout << "We will cut time range in " << npts+1 << " intervals of " << ndays << " days " << std::endl;
    
  //
  // Maintenant on peut definir les graphiques qui vons contenir les données
  //
    
  std::vector<TH2F* > SC_vs_time;
  std::vector<TH2F* > SC_vs_time_raw;
  std::vector<TH2F* > Ratio_vs_time;
  std::vector<TH2F* > Ratio_vs_time_raw;
    
  std::vector<std::vector<std::pair<double, double> > > totaldata;
  std::vector<std::pair<double, double> > SCdate;
  SCdate.clear();
  std::pair<double, double> couple;

  for (int i=0;i<50;++i)
  {
      TH2F *aplotxy  = new TH2F("example plot 2Da","a",  2*npts_r,tmin,tmaxreal, 100,0.,1.1*valmax);
      TH2F *raw      = new TH2F("example plot 2Db","b", 2*npts_r,tmin,tmaxreal, 200,0.,1.1*rawvalmax);
      TH2F *ratio    = new TH2F("example plot 1Dc","c", 2*npts_r,tmin,tmaxreal,400,0.,40.);
      TH2F *ratior   = new TH2F("example plot 1Dd","d", 2*npts_r,tmin,tmaxreal,250,0.,25.);
      
      SC_vs_time.push_back(aplotxy);
      SC_vs_time_raw.push_back(raw);
      Ratio_vs_time.push_back(ratio);
      Ratio_vs_time_raw.push_back(ratior);
      totaldata.push_back(SCdate);
  }
   
  TH2F *Ratio_vs_age    = new TH2F("example plot 2Dd","d", ndays*(npts+1),tmin,tmax, 4,20.,100);
  TH2F *Ratioraw_vs_age = new TH2F("example plot 2De","e", ndays*(npts+1),tmin,tmax, 4,20.,100);
    
  //
  // Deuxième boucle sur les données, on collecte et on met en forme
  //
    
  for (int i=0;i<nentries;++i)
  {
    newtree->GetEntry(i);
    if (vac>6) continue;
    if (pop<=minbycat) continue;
    if (age==2) continue;  // Catégorie pas étudiée
    
    couple.first=date;
    couple.second=val+1e12*pop;

    totaldata[int(10*age+vac)].push_back(couple);
  
    // Ici on additione toutes les données de vaccinés complets
    if (vac<6)
    {
        totaldata[int(10*age+7)].push_back(couple);
        totaldata[8].push_back(couple); // Toutes les données vax
    }
    else
    {
        totaldata[9].push_back(couple); // Toutes les données non-vax
    }
      
  }

  //
  // On a les données pour chaque catégories dans des vecteurs
  // on les classes par date croissante et on fait les plots
  //
    
  for (int i=0;i<50;++i)
  {
      sort(totaldata[i].begin(), totaldata[i].end(), myComparison);
      doPlot(&totaldata[i],SC_vs_time[i],SC_vs_time_raw[i],tmin,tmax,ndays);
  }

  int age_r[5] = {50,90,10,30,70};
    
  for (int i=0;i<5;++i)
      doRatioPlot(&totaldata[10*i+7],&totaldata[10*i+6],Ratio_vs_time[i],Ratio_vs_time_raw[i],age_r[i],Ratio_vs_age,Ratioraw_vs_age,tmin,tmax,ndays);

   doRatioPlot(&totaldata[8],&totaldata[9],Ratio_vs_time[8],Ratio_vs_time_raw[8],-10,Ratio_vs_age,Ratioraw_vs_age,tmin,tmax,ndays);
  //
  // Les plots sont finis on fait la cosmétique
  //

  TCanvas *c1 = new TCanvas("c1","Augmentation du risque d'admission (nonvax/vax)",451,208,1208,604);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetBorderSize(2);
  c1->SetLeftMargin(0.08);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderMode(0);

  TLegend *leg = new TLegend(0.1,0.65,0.3,0.9);

  leg->SetTextSize(0.035);
  leg->SetHeader("Tranches d'#hat{a}ge","C");
  leg->SetFillColor(0);
  leg->AddEntry(Ratio_vs_time[3],"[20,39]","p");
  leg->AddEntry(Ratio_vs_time[0],"[40,59]","p");
  leg->AddEntry(Ratio_vs_time[4],"[60,79]","p");
  leg->AddEntry(Ratio_vs_time[1],"[80+]","p");

  Ratio_vs_time[3]->SetMarkerColor(38);
  Ratio_vs_time[0]->SetMarkerColor(28);
  Ratio_vs_time[4]->SetMarkerColor(1);
  Ratio_vs_time[1]->SetMarkerColor(2);
    
  Ratio_vs_time[0]->GetXaxis()->SetTimeDisplay(1);
  Ratio_vs_time[0]->GetXaxis()->SetTimeFormat("%d/%m/%Y");
  Ratio_vs_time[0]->GetXaxis()->SetTimeOffset(0,"gmt");
  Ratio_vs_time[0]->GetXaxis()->SetLabelOffset(0.015);
  Ratio_vs_time[0]->GetXaxis()->SetLabelSize(0.03);
  Ratio_vs_time[0]->GetYaxis()->SetLabelSize(0.03);
  Ratio_vs_time[0]->GetXaxis()->SetNdivisions(-507);
  Ratio_vs_time[0]->GetYaxis()->SetTitle("Augmentation du risque d'admission (nonvax/vax)");
  Ratio_vs_time[0]->SetTitleSize(0.035);
  Ratio_vs_time[0]->SetMarkerStyle(20);
  Ratio_vs_time[0]->Draw();
  
  for (int i=1;i<5;++i)
  {
    if (i==2) continue;
    Ratio_vs_time[i]->SetMarkerStyle(20);
    Ratio_vs_time[i]->Draw("same");
  }
  leg->Draw();
    
  c1->Update();
    
    
    TCanvas *c11 = new TCanvas("c11","Augmentation du risque d'admissions (nonvax/vax)",451,208,1208,604);
    c11->SetFillColor(0);
    c11->SetBorderMode(0);
    c11->SetGridx();
    c11->SetGridy();
    c11->SetBorderSize(2);
    c11->SetLeftMargin(0.08);
    c11->SetFrameBorderMode(0);
    c11->SetFrameBorderMode(0);
            
    Ratio_vs_time[8]->GetXaxis()->SetTimeDisplay(1);
    Ratio_vs_time[8]->GetXaxis()->SetTimeFormat("%d/%m/%Y");
    Ratio_vs_time[8]->GetXaxis()->SetTimeOffset(0,"gmt");
    Ratio_vs_time[8]->GetXaxis()->SetLabelOffset(0.015);
    Ratio_vs_time[8]->GetXaxis()->SetLabelSize(0.03);
    Ratio_vs_time[8]->GetYaxis()->SetLabelSize(0.03);
    Ratio_vs_time[8]->GetXaxis()->SetNdivisions(-507);
    Ratio_vs_time[8]->GetYaxis()->SetTitle("Augmentation du risque d'admissions (nonvax/vax)");
    Ratio_vs_time[8]->SetTitleSize(0.035);
    Ratio_vs_time[8]->SetMarkerStyle(20);
    Ratio_vs_time[8]->Draw();
      
    c11->Update();

    TCanvas *c11b = new TCanvas("c11b","Rapport du nombre d'admissions (nonvax/vax)",451,208,1208,604);
    c11b->SetFillColor(0);
    c11b->SetBorderMode(0);
    c11b->SetGridx();
    c11b->SetGridy();
    c11b->SetBorderSize(2);
    c11b->SetLeftMargin(0.08);
    c11b->SetFrameBorderMode(0);
    c11b->SetFrameBorderMode(0);
            
    Ratio_vs_time_raw[8]->GetXaxis()->SetTimeDisplay(1);
    Ratio_vs_time_raw[8]->GetXaxis()->SetTimeFormat("%d/%m/%Y");
    Ratio_vs_time_raw[8]->GetXaxis()->SetTimeOffset(0,"gmt");
    Ratio_vs_time_raw[8]->GetXaxis()->SetLabelOffset(0.015);
    Ratio_vs_time_raw[8]->GetXaxis()->SetLabelSize(0.03);
    Ratio_vs_time_raw[8]->GetYaxis()->SetLabelSize(0.03);
    Ratio_vs_time_raw[8]->GetXaxis()->SetNdivisions(-507);
    Ratio_vs_time_raw[8]->GetYaxis()->SetTitle("Rapport du nombre d'admissions (nonvax/vax)");
    Ratio_vs_time_raw[8]->SetTitleSize(0.035);
    Ratio_vs_time_raw[8]->SetMarkerStyle(20);
    Ratio_vs_time_raw[8]->Draw();
      
    c11b->Update();
    
    
  TCanvas *c3 = new TCanvas("c3","Non-vaccines, donnees par millions d'habitants ",451,208,1208,604);
  c3->SetFillColor(0);
  c3->SetBorderMode(0);
  c3->SetGridx();
  c3->SetGridy();
  c3->SetBorderSize(2);
  c3->SetLeftMargin(0.08);
  c3->SetFrameBorderMode(0);
  c3->SetFrameBorderMode(0);
      
  TLegend *leg3 = new TLegend(0.1,0.65,0.3,0.9);

  leg3->SetTextSize(0.035);
  leg3->SetHeader("Tranches d'#hat{a}ge","C");
  leg3->SetFillColor(0);
  leg3->AddEntry(SC_vs_time[36],"[20,39]","p");
  leg3->AddEntry(SC_vs_time[6],"[40,59]","p");
  leg3->AddEntry(SC_vs_time[46],"[60,79]","p");
  leg3->AddEntry(SC_vs_time[16],"[80+]","p");

  SC_vs_time[36]->SetMarkerColor(38);
  SC_vs_time[6]->SetMarkerColor(28);
  SC_vs_time[46]->SetMarkerColor(1);
  SC_vs_time[16]->SetMarkerColor(2);
      
  SC_vs_time[6]->GetXaxis()->SetTimeDisplay(1);
  SC_vs_time[6]->GetXaxis()->SetTimeFormat("%d/%m/%Y");
  SC_vs_time[6]->GetXaxis()->SetTimeOffset(0,"gmt");
  SC_vs_time[6]->GetXaxis()->SetLabelOffset(0.015);
  SC_vs_time[6]->GetXaxis()->SetLabelSize(0.03);
  SC_vs_time[6]->GetYaxis()->SetLabelSize(0.03);
  SC_vs_time[6]->GetXaxis()->SetNdivisions(-507);
  SC_vs_time[6]->GetYaxis()->SetTitle("Nombre moyen d'admissions par jour (par millions d'habitants)");
  SC_vs_time[6]->SetTitleSize(0.035);
  SC_vs_time[6]->SetMarkerStyle(20);
  SC_vs_time[6]->Draw();
    
  for (int i=1;i<5;++i)
  {
    if (i==2) continue;
    SC_vs_time[10*i+6]->SetMarkerStyle(20);
    SC_vs_time[10*i+6]->Draw("same");
  }
  leg3->Draw();
    
  c3->Update();
 
    
  TCanvas *c4 = new TCanvas("c4","Non vaccines, données brutes",451,208,1208,604);
  c4->SetFillColor(0);
  c4->SetBorderMode(0);
    c4->SetGridx();
    c4->SetGridy();
    c4->SetBorderSize(2);
    c4->SetLeftMargin(0.08);
    c4->SetFrameBorderMode(0);
    c4->SetFrameBorderMode(0);
      
    TLegend *leg4 = new TLegend(0.1,0.65,0.3,0.9);

    leg4->SetTextSize(0.035);
    leg4->SetHeader("Tranches d'#hat{a}ge","C");
    leg4->SetFillColor(0);
    leg4->AddEntry(SC_vs_time[36],"[20,39]","p");
    leg4->AddEntry(SC_vs_time[6],"[40,59]","p");
    leg4->AddEntry(SC_vs_time[46],"[60,79]","p");
    leg4->AddEntry(SC_vs_time[16],"[80+]","p");

      SC_vs_time_raw[36]->SetMarkerColor(38);
      SC_vs_time_raw[6]->SetMarkerColor(28);
      SC_vs_time_raw[46]->SetMarkerColor(1);
      SC_vs_time_raw[16]->SetMarkerColor(2);
      
    SC_vs_time_raw[6]->GetXaxis()->SetTimeDisplay(1);
    SC_vs_time_raw[6]->GetXaxis()->SetTimeFormat("%d/%m/%Y");
    SC_vs_time_raw[6]->GetXaxis()->SetTimeOffset(0,"gmt");
    SC_vs_time_raw[6]->GetXaxis()->SetLabelOffset(0.015);
    SC_vs_time_raw[6]->GetXaxis()->SetLabelSize(0.03);
    SC_vs_time_raw[6]->GetYaxis()->SetLabelSize(0.03);
    SC_vs_time_raw[6]->GetXaxis()->SetNdivisions(-507);
    SC_vs_time_raw[6]->GetYaxis()->SetTitle("Nombre moyen d'admissions par jour");
    SC_vs_time_raw[6]->SetTitleSize(0.035);
    SC_vs_time_raw[6]->SetMarkerStyle(20);
    SC_vs_time_raw[6]->Draw();
    
    for (int i=1;i<5;++i)
    {
      if (i==2) continue;
      SC_vs_time_raw[10*i+6]->SetMarkerStyle(20);
      SC_vs_time_raw[10*i+6]->Draw("same");
    }
    leg4->Draw();
      
    c4->Update();

    TCanvas *c7 = new TCanvas("c7","Schema vaccinal complet, donnees brutes",451,208,1208,604);
    c7->SetFillColor(0);
    c7->SetBorderMode(0);
    c7->SetGridx();
    c7->SetGridy();
    
    c7->SetBorderSize(2);
    c7->SetLeftMargin(0.08);
    c7->SetFrameBorderMode(0);
    c7->SetFrameBorderMode(0);
        
    TLegend *leg7 = new TLegend(0.1,0.65,0.3,0.9);
    leg7->SetTextSize(0.035);
    leg7->SetHeader("Tranches d'#hat{a}ge","C");
    leg7->SetFillColor(0);
    leg7->AddEntry(SC_vs_time_raw[37],"[20,39]","p");
    leg7->AddEntry(SC_vs_time_raw[7],"[40,59]","p");
    leg7->AddEntry(SC_vs_time_raw[47],"[60,79]","p");
    leg7->AddEntry(SC_vs_time_raw[17],"[80+]","p");

    SC_vs_time_raw[37]->SetMarkerColor(38);
    SC_vs_time_raw[7]->SetMarkerColor(28);
    SC_vs_time_raw[47]->SetMarkerColor(1);
    SC_vs_time_raw[17]->SetMarkerColor(2);
        
    SC_vs_time_raw[7]->GetXaxis()->SetTimeDisplay(1);
    SC_vs_time_raw[7]->GetXaxis()->SetTimeFormat("%d/%m/%Y");
    SC_vs_time_raw[7]->GetXaxis()->SetTimeOffset(0,"gmt");
    SC_vs_time_raw[7]->GetXaxis()->SetLabelOffset(0.015);
    SC_vs_time_raw[7]->GetXaxis()->SetLabelSize(0.03);
    SC_vs_time_raw[7]->GetYaxis()->SetLabelSize(0.03);
    SC_vs_time_raw[7]->GetXaxis()->SetNdivisions(-507);
    SC_vs_time_raw[7]->GetYaxis()->SetTitle("Nombre moyen d'admissions par jour");
    SC_vs_time_raw[7]->SetTitleSize(0.035);
    SC_vs_time_raw[7]->SetMarkerStyle(20);
    SC_vs_time_raw[7]->Draw();
      
    for (int i=1;i<5;++i)
    {
     if (i==2) continue;
     SC_vs_time_raw[10*i+7]->SetMarkerStyle(20);
     SC_vs_time_raw[10*i+7]->Draw("same");
    }
    
    leg7->Draw();
    c7->Update();
    
    
    TCanvas *c8 = new TCanvas("c8","Schema vaccinal complet, donnees par millions d'habitants",451,208,1208,604);
      c8->SetFillColor(0);
      c8->SetBorderMode(0);
      c8->SetGridx();
      c8->SetGridy();
      c8->SetBorderSize(2);
      c8->SetLeftMargin(0.08);
      c8->SetFrameBorderMode(0);
      c8->SetFrameBorderMode(0);
        
      TLegend *leg8 = new TLegend(0.1,0.65,0.3,0.9);

      leg8->SetTextSize(0.035);
      leg8->SetHeader("Tranches d'#hat{a}ge","C");
      leg8->SetFillColor(0);
      leg8->AddEntry(SC_vs_time[37],"[20,39]","p");
      leg8->AddEntry(SC_vs_time[7],"[40,59]","p");
      leg8->AddEntry(SC_vs_time[47],"[60,79]","p");
      leg8->AddEntry(SC_vs_time[17],"[80+]","p");

      SC_vs_time[37]->SetMarkerColor(38);
      SC_vs_time[7]->SetMarkerColor(28);
      SC_vs_time[47]->SetMarkerColor(1);
      SC_vs_time[17]->SetMarkerColor(2);
        
      SC_vs_time[7]->GetXaxis()->SetTimeDisplay(1);
      SC_vs_time[7]->GetXaxis()->SetTimeFormat("%d/%m/%Y");
      SC_vs_time[7]->GetXaxis()->SetTimeOffset(0,"gmt");
      SC_vs_time[7]->GetXaxis()->SetLabelOffset(0.015);
      SC_vs_time[7]->GetXaxis()->SetLabelSize(0.03);
      SC_vs_time[7]->GetYaxis()->SetLabelSize(0.03);
      SC_vs_time[7]->GetXaxis()->SetNdivisions(-507);
      SC_vs_time[7]->GetYaxis()->SetTitle("Nombre moyen d'admissions par jour (par millions d'habitants)");
      SC_vs_time[7]->SetTitleSize(0.035);
      SC_vs_time[7]->SetMarkerStyle(20);
      SC_vs_time[7]->Draw();
      
      for (int i=1;i<5;++i)
      {
        if (i==2) continue;
        SC_vs_time[10*i+7]->SetMarkerStyle(20);
        SC_vs_time[10*i+7]->Draw("same");
      }
      leg8->Draw();
        
      c8->Update();

     
    gStyle->SetPalette(55);
    gStyle->SetNumberContours(256);
    TCanvas *c9 = new TCanvas("c9","Rapports normalisé par millions d'habitants non-vax/vax en fonction de l'age pour la variable",451,208,1208,604);
    c9->SetFillColor(0);
    c9->SetBorderMode(0);
    c9->SetBorderSize(2);
    //c9->Range(1.620928e+09,9.999999,1.639473e+09,110);
    c9->SetLogz(1);
    c9->SetLeftMargin(0.08);
    c9->SetFrameBorderMode(0);
    c9->SetFrameBorderMode(0);

    Ratio_vs_age->GetXaxis()->SetTimeDisplay(1);
    Ratio_vs_age->GetXaxis()->SetTimeFormat("%d/%m/%Y");
    Ratio_vs_age->GetXaxis()->SetTimeOffset(0,"gmt");
    Ratio_vs_age->GetXaxis()->SetLabelOffset(0.015);
    Ratio_vs_age->GetXaxis()->SetLabelSize(0.03);
    Ratio_vs_age->GetYaxis()->SetLabelSize(0.03);
    Ratio_vs_age->GetXaxis()->SetNdivisions(-507);
    Ratio_vs_age->GetYaxis()->SetTitle("Age des patients");
    Ratio_vs_age->SetTitleSize(0.035);
    Ratio_vs_age->SetMarkerStyle(20);
    Ratio_vs_age->Draw("colz");
    
    c9->Update();
    
    TCanvas *c10 = new TCanvas("c10","Rapports bruts non-vax/vax en fonction de l'age pour la variable",451,208,1208,604);
    c10->SetFillColor(0);
    c10->SetBorderMode(0);
    c10->SetBorderSize(2);
    c10->SetLogz(1);
    c10->SetLeftMargin(0.08);
    c10->SetFrameBorderMode(0);
    c10->SetFrameBorderMode(0);
    
    Ratioraw_vs_age->GetXaxis()->SetTimeDisplay(1);
    Ratioraw_vs_age->GetXaxis()->SetTimeFormat("%d/%m/%Y");
    Ratioraw_vs_age->GetXaxis()->SetTimeOffset(0,"gmt");
    Ratioraw_vs_age->GetXaxis()->SetLabelOffset(0.015);
    Ratioraw_vs_age->GetXaxis()->SetLabelSize(0.03);
    Ratioraw_vs_age->GetYaxis()->SetLabelSize(0.03);
    Ratioraw_vs_age->GetXaxis()->SetNdivisions(-507);
    Ratioraw_vs_age->GetYaxis()->SetTitle("Age des patients");
    Ratioraw_vs_age->SetTitleSize(0.035);
    Ratioraw_vs_age->SetMarkerStyle(20);
    Ratioraw_vs_age->Draw("colz");
    
    c10->Update();
}



