{
	gROOT->Reset();
	gStyle->SetPalette(1);
	gStyle->SetOptFit(11111);
	// gStyle->SetPadRightMargin(0.32);

	//Draw a simple graph
	// To see the output of this macro, click begin_html <a href="gif/graph.gif">here</a>. end_html
	//Author: Rene Brun

	TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);

	//    c1->SetFillColor(42);
	//    c1->SetGrid();

	const Int_t n = 9;
	const Int_t n2=3;
	Double_t E[n]={2.0,2.6,3.2,3.8,4.4,5,5.4,6,6.8};
	Double_t E_error[n]={0,0,0,0,0,0,0,0,0};   
	Double_t E2[n2]={5,6,7};
	Double_t E2_error[n2]={0,0,0};   
	Double_t Mean[n]={1.97827e+03,2.57895e+03,3.18330e+03,3.76794e+03,4.39214e+03,4.97318e+03,5.37839e+03,5.97958e+03,6.78001e+03};
	Double_t Mean_error[n]={1.17850e+00,1.38842e+00,1.60310e+00,1.67957e+00,1.74786e+00,1.88027e+00,2.00896e+00,2.08335e+00,2.21579e+00};
	Double_t Sigma[n]={7.41022e+01,8.54021e+01,9.07927e+01,9.84583e+01,1.02737e+02,1.09864e+02,1.15787e+02,1.20565e+02,1.25743e+02};
	Double_t Sigma_error[n]={8.92351e-01,1.13155e+00,1.41997e+00,1.26777e+00,1.25459e+00,1.46807e+00,1.53518e+00,1.68004e+00,1.63854e+00};
	// all angle no field 25degree Calibrate
	Double_t Mean1[n]={2.00292e+03,2.61406e+03,3.19296e+03,3.79870e+03,4.40946e+03,5.01457e+03,5.42501e+03,6.03050e+03,6.83945e+03};
	Double_t Mean1_error[n]={1.81052e+00,1.46137e+00,1.57466e+00,1.68476e+00,1.77712e+00,1.91881e+00,2.04070e+00,2.12692e+00,2.37970e+00};
	Double_t Sigma1[n]={8.01272e+01,9.16308e+01,9.96573e+01,1.05586e+02,1.13418e+02,1.18919e+02,1.22648e+02,1.27819e+02,1.36128e+02};
	Double_t Sigma1_error[n]={1.44504e+00,1.07504e+00,1.22052e+00,1.20927e+00,1.37034e+00,1.47494e+00,1.54418e+00,1.65008e+00,1.85036e+00};

	//straight hit

	//nofield_all angle_Ecaliberate

	Double_t Mean2[n]={1.99812e+03,2.59874e+03,3.20220e+03,3.80095e+03,4.40239e+03,5.00537e+03,5.40310e+03,6.00540e+03,6.80430e+03};
	Double_t Mean2_error[n]={1.47823e+00,1.65066e+00,1.84160e+00,1.95710e+00,2.11863e+00,2.32843e+00,2.42389e+00,2.46719e+00,3.03493e+00};
	Double_t Sigma2[n]={9.00860e+01,9.92363e+01,1.07110e+02,1.11866e+02,1.19338e+02,1.29375e+02,1.32591e+02,1.39075e+02,1.43973e+02};
	Double_t Sigma2_error[n]={1.12043e+00,1.21856e+00,1.40164e+00,1.38704e+00,1.67091e+00,1.85352e+00,1.82587e+00,1.93933e+00,2.31251e+00};

	// preshower+shower theta=35 calibrate
	Double_t Mean3[n]={1.95245e+03,2.47693e+03,3.07127e+03,3.66925e+03,4.26781e+03,4.86791e+03,5.29002e+03,5.86752e+03,6.69130e+03};
	Double_t Mean3_error[n]={4.53010e-01,6.05287e-01,8.22281e-01,9.93725e-01,1.23506e+00,7.57264e-01,1.54250e+00,1.66642e+00,1.60568e+00};
	Double_t Sigma3[n]={9.80468e+01,1.07549e+02,1.16086e+02,1.23163e+02,1.32411e+02,1.38464e+02,1.39807e+02,1.46634e+02,1.55908e+02};
	Double_t Sigma3_error[n]={3.56171e-01,4.43996e-01,7.48182e-01,9.20442e-01,1.28177e+00,8.38823e-01,1.58278e+00,1.55366e+00,1.20708e+00};
	//6+1cluster test
	Double_t Mean4[n2]={1.10503e+03,6.10199e+02,5.06705e+02};
	Double_t Mean4_error[n2]={9.36402e-01,1.38512e+01,1.53274e+00};
	Double_t Sigma4[n2]={4.49333e+01,3.52025e+02,3.53324e+02};
	Double_t Sigma4_error[n2]={7.29055e-01,1.32737e+01,9.80118e-01};

	Double_t Mean5[n2]={1.07615e+03,4.47071e+02,3.89256e+02};
	Double_t Mean5_error[n2]={1.73061e+00,2.78417e+01,2.05917e+00};
	Double_t Sigma5[n2]={3.99675e+01,4.18375e+02,3.66666e+02};
	Double_t Sigma5_error[n2]={1.33466e+00,2.69196e+01,1.20556e+00};

	Double_t res[n];
	Double_t res_error[n];   
	Double_t factor[n];
	Double_t res1[n];
	Double_t res1_error[n]; 
	//factor[n]=Mean[n]/E[n];
	Double_t res2[n];
	Double_t res2_error[n]; 
	Double_t res3[n];
	Double_t res3_error[n]; 
	Double_t res4[n];
	Double_t res4_error[n]; 
	Double_t res5[n];
	Double_t res5_error[n]; 
	for (Int_t i=0;i<n2;i++) {
		res4[i]=Mean4[i]/1e3/E2[i];
		res4_error[i]=Mean4_error[i]/1e3/E2[i];
		res5[i]=Mean5[i]/1e3/E2[i];
		res5_error[i]=Mean5_error[i]/1e3/E2[i];
	}

	for (Int_t i=0;i<n;i++) {
		res[i]=Sigma[i]/Mean[i];
		res_error[i]=Sigma_error[i]/Mean[i]; 
		res1[i]=Sigma1[i]/Mean1[i];;
		res1_error[i]=Sigma1_error[i]/Mean1[i];
		res2[i]=Sigma2[i]/Mean2[i];;
		res2_error[i]=Sigma2_error[i]/Mean2[i];
		res3[i]=Sigma3[i]/Mean3[i];;
		res3_error[i]=Sigma3_error[i]/Mean3[i];
	}
	gr = new TGraphErrors(n,E,res,E_error,res_error);
	gr1 = new TGraphErrors(n,E,res1,E_error,res1_error);
	gr2 = new TGraphErrors(n,E,res2,E_error,res2_error);
	gr3 = new TGraphErrors(n,E,res3,E_error,res3_error);

	gr4 = new TGraphErrors(n2,E2,res4,E2_error,res4_error);
	gr5 = new TGraphErrors(n2,E2,res5,E2_error,res5_error);
	gr->Print();
	gr->SetMarkerSize(1.5);
	gr->SetMarkerStyle(22);
	gr->SetMarkerColor(2);
	gr1->SetTitle("EC calibrated energy(shower+preshower) / E_toal");
	gr1->GetXaxis()->SetTitle("E (GeV)");
	gr1->GetYaxis()->SetTitle("Edep/E");

	//TF1 *fun = new TF1("fun","sqrt(pow([0]/sqrt(x),2)+pow([1],2))",1,7);
	TF1 *fun = new TF1("fun","sqrt(pow([0]/sqrt(x),2)+pow([1],2)+pow([2]/x,2))",1,7);
	fun->SetLineStyle(9);
	fun->SetLineColor(4);
	TF1 *fun2 = new TF1("fun2","sqrt(pow([0]/sqrt(x),2)+pow([1],2)+pow([2]/x,2))",1,7);
	fun2->SetLineStyle(9);
	fun2->SetLineColor(6);
	TF1 *fun3 = new TF1("fun3","sqrt(pow([0]/sqrt(x),2)+pow([1],2)+pow([2]/x,2))",1,7);
	fun3->SetLineStyle(9);
	fun3->SetLineColor(2);
	gr1->SetMarkerSize(1.5);
	gr1->SetMarkerStyle(22);
	gr1->SetMarkerColor(6);
	gr1->SetMaximum(0.07);
	gr1->Fit("fun2","R");
	gr1->Draw("Ap"); 
	gr2->SetMarkerSize(1.);
	gr2->SetMarkerColor(1);
	gr2->SetMarkerStyle(21);
	gr2->Fit("fun","R");
	gr2->Draw("same p"); 
	gr3->SetMarkerSize(1.);   
	gr3->SetMarkerColor(4);
	gr3->SetMarkerStyle(22);
	gr3->Draw("same p");
	gr3->Fit("fun","R+");
	gr4->SetMarkerSize(1.);
	gr4->SetMarkerColor(4);
	gr4->SetMarkerStyle(22);
	gr5->SetMarkerSize(1.);
	gr5->SetMarkerColor(4);
	gr5->SetMarkerStyle(22); 

	TLegend *leg2 = new TLegend(0.11,0.71,0.30,0.89,NULL,"brNDC");//0.6,0.78,0.9,0.90
	leg2->AddEntry(gr1,"nofield 25#circ PVDIS FAEC vertex(0,0,10cm)","p");
	leg2->AddEntry(gr3,"nofield 35#circ PVDIS FAEC vertex(0,0,10cm)","p");
	leg2->AddEntry(gr2,"nofield [22,35]#circ PVDIS FAEC vertex(0,0,10cm)","p");
	leg2->SetBorderSize(0);
	leg2->SetTextFont(62);
	leg2->SetTextSize(0.03);
	leg2->SetFillColor(0);
	leg2->SetFillStyle(1001);
	leg2->Draw();

}
