#include <vector>
#include <class_B_field2.cpp>
#include <fstream>
#include <iostream>
#include <boost/lexical_cast.hpp>
#include <string>
#include <iomanip>

int bfield_out(int fname_num, int ftype){
 //Conditional for field type: 11=output regular field, 7=output KRF
 if(ftype==7){
    string field_outfile = "krf_field_out_";
    string fname_str = boost::lexical_cast<string>(fname_num);
    field_outfile.append(fname_str);
    field_outfile.append(".txt");
    ofstream bfout;
    bfout.open(field_outfile.c_str());
    for (uint i=0;i<brfx.size();i++){
    bfout << std::fixed;    
    bfout << setprecision(15) << brfx[i] << '\t' << brfy[i] << '\t'
    << brfz[i] << setprecision(3) << '\t' << brxc[i] << '\t' 
    << bryc[i] << '\t' << brzc[i] << std::endl;
    }
    bfout.close();
 }
 if(ftype==11){
    string field_outfile = "reg_field_out_";
    string fname_str = boost::lexical_cast<string>(fname_num);
    field_outfile.append(fname_str);
    field_outfile.append(".txt");
    ofstream bfout;
    bfout.open(field_outfile.c_str());
    for (uint i=0;i<bregfx.size();i++){
    bfout << std::fixed;    
    bfout << setprecision(15) << bregfx[i] << '\t' << bregfy[i] << '\t'
    << bregfz[i] << setprecision(3) << '\t' << bregxc[i] << '\t' 
    << bregyc[i] << '\t' << bregzc[i] << std::endl;
    }
    bfout.close();
 }
return 0;
}
