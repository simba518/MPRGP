#ifndef _PYTHONFFT_H_
#define _PYTHONFFT_H_

namespace UTILITY{

  class PythonFFT{
  
  public:
	template <class VECTOR>
	void add(const string name,const VECTOR& s,const double h){
	  stringstream data;
	  data << name<<"=[";
	  for (int i = 0; i < (int)s.size(); ++i){
		data << std::setprecision(14) << s[i];
		if((int)s.size()-1!=i) data<<", ";
	  }
	  data << "]\n";
	  data << name<<"_h="<<h<<"\n\n";
	  signalData += data.str();
	  if(signalList.size()>0){
		signalList = signalList+","+name;
		stepList = stepList+","+name+"_h";
	  }else{
		stepList = name+"_h";
		signalList = name;
	  }
	}
	void write(const string fname,const string saveFigTo=""){
	  OUTFILE(file,string(fname+".py").c_str());
	  file << head();
	  file << signalData;
	  file << "\nsignalList=[" << signalList <<"]\n";
	  file << "\nstepList=[" << stepList <<"]\n";
	  file << end(saveFigTo);
	}
	void clear(){
	  signalData = "";
	  signalList = "";
	  stepList = "";
	}
  
  protected:
	static string head(){
	  string h = "import numpy as np\n";
	  h += "import scipy\n";
	  h += "import scipy.fftpack\n";
	  h += "import matplotlib\n";
	  h += "import matplotlib.pyplot as plt\n";
	  h += "import sys\n\n";
	  return h;
	}
	static string end(const string saveFigTo){
	  string e = "\nfor i in range(0,len(signalList)):\n";
	  e += "\tsignal = signalList[i]\n";
	  e += "\tfourier = np.fft.fft(signal)\n";
	  e += "\tfreq = np.fft.fftfreq(len(signal),stepList[i])\n";
	  e += "\tplt.stem(freq,np.absolute(fourier),\"-\")\n";
	  e += "\tplt.grid(True)\n";
	  if (saveFigTo.size()>0){
		e += string("\tplt.savefig(\"")+saveFigTo+"_\"+str(i)+\".png\")\n";
	  }else{
		e += "\tplt.show()\n";
	  }
	  e += "\tplt.clf()";
	  return e;
	}
  
  private:
	string signalData;
	string signalList;
	string stepList;
  };
}

#endif /* _PYTHONFFT_H_ */

