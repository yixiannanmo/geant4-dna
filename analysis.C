
{
  gROOT->Reset();//重置堆栈ROOT说明书P151
 
  gStyle->SetOptStat("em");//选择显示信息的类型，见root说明书P77，e表示entries，m表示mean
  TCanvas *c1;
  TPad *pad1, *pad2, *pad3;//利用TPad画几个子pad，并且画出来P178
  c1 = new TCanvas("c1","PDB DNA outputs",200,10,700,780);
  c1->SetFillColor(0);

  pad1 = new TPad("pad1","pad1",0.02,0.52,0.98,0.98,21);//都类似于Tcanvas（是Tpad的子类）
  pad2 = new TPad("pad2","pad2",0.02,0.02,0.48,0.48,21);
  pad3 = new TPad("pad3","pad3",0.52,0.02,0.98,0.48,21);

  pad1->SetFillColor(0);
  pad1->Draw();
  pad2->SetFillColor(0);
  pad2->Draw();
  pad3->SetFillColor(0);
  pad3->Draw();

//引入G4得到的数据
  TFile f("pdb4dna_output.root");//利用Tfile读取数据后面还可以加参数（“recreate”）


  TH1D* hist1 = (TH1D*)f.Get("1");//TFile::Get("Tkey关键词来调用每个pad")
  pad1->cd();//Tpad：：cd（）激活画布
  hist1->Draw("HIST");

  TH1D* hist2 = (TH1D*)f.Get("2");
  pad2->cd();
  hist2->Draw("HIST");

  TH1D* hist3 = (TH1D*)f.Get("3");
  pad3->cd();
  hist3->Draw("HIST");

  c1->Modified();//画布改变，重新调整P188
  c1->Update();



  double* pdbStats=new double[4];

  hist1->GetStats(pdbStats);
  cout << "-> Edep in the target : " << pdbStats[2]/1E6 << " MeV" << endl;

  hist2->GetStats(pdbStats);
  cout << "-> Number of SSB : " << pdbStats[2] << endl;

  hist3->GetStats(pdbStats);
  cout << "-> Number of DSB : " << pdbStats[2] << endl;
}
