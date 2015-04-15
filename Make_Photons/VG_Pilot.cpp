// version: 6.1.0    date: 2014-05-06


/*    VG_Pilot.cpp class
         Charlie Duke
        August 28, 2010
        Can extract matlab formated tokens

*/
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <map>

using namespace std;

// needed for extern declaration for oLog output stream
//#include "VGA_Definition.h"

#include "VG_Pilot.h"

#define DEBUG(x) *oLog << #x << " = " << x << endl

//******************************************
VG_Pilot::VG_Pilot() {

}

//*************** end of VG_Pilot ****************
VG_Pilot::VG_Pilot (const string filename) {

  // the addPilotFile method loads the maps
  addPilotFile(filename);
}

//************** end of VG_Pilot ****************

bool VG_Pilot::addPilotFile(const string &pilotfile) {
  char chB[] = " ";  // blank character
  char chA[] = "*";  // asterisk character

  // open pilot file for input

  int debug = 1;
  if (debug >0) {
    *oLog << " -- VG_Pilot::addPilotFile " << endl;
    *oLog << "       adding pilotfile: " << pilotfile << endl;
  }

  ifstream inFile(pilotfile.c_str(),ios::in);
  if (! inFile) {
    cerr << "error: not able to open pilot file: "
         << pilotfile << endl;
    exit(0);
  }
  string pilotline;
  while( getline(inFile,pilotline,'\n')) {

    // skip blank lines
    if (pilotline.size() == 0) continue;

    // see if first character is an "*"
    std::string::size_type idx = 0;
    idx = pilotline.find_first_not_of(chB);

    // skip line containing only whitespace
    if (idx == string::npos) continue;

    // see if we have a leading asterisk
    if (pilotline.compare(idx,1,chA,1)==0) {
      // we have an asterisk, remove it and tokenize the line
      
      pilotline.erase(idx,1);  

      tokenizer(pilotline, tokens);
      //cerr << "tokenssize: " << tokens.size() << endl;

      if (tokens.size() == 0) continue;
       
      flag = tokens[0];
      // erase tokens[0] and store the rest in the map
      tokens.erase(tokens.begin());
      entry p1(flag,tokens);
      m.insert(p1);
      
    }
    
  }

  return true;
}

//********************* end of addPilotFile ******************8

int VG_Pilot::set_flag(const string flagname) {
  // set variables to access the first flag line
  flag = flagname;
  line_count = 0;
  count = m.count(flag);

  if (count) {
    lower = m.lower_bound(flag);
    upper = m.upper_bound(flag);
  }
  iter = lower;
  return count;
}
//**************** end of set_flag ***************************

int VG_Pilot::get_line_vector(vector<string> &line_vector) {

  line_vector.clear();

  vector<string>::iterator viter;

  if (count>0) {
    line_vector = (*iter).second;
    iter++;
  }
  //else {  
  //line_vector.clear();
  //}
  count--;  // decrement count
  return count;
}
//******************** end of get_line_vector

void VG_Pilot::tokenizer(string& str, vector<string>& tokens) {

  tokens.clear();
  string strc = str;
  // string can contain matlab entries, e.g. [1 3 5], so simple
  // tokenizer doesn't work.

  char chB[] = " ";
  char chL[] = "[";
  char chR[] = "]";
  string::size_type lgt = 1;

  string::size_type idx = 0;

  // find first non-blank char starting at idx location
  while ( (idx = strc.find_first_not_of(chB,idx,lgt) ) !=string::npos ) {

    // set start position of token
    string::size_type start = idx;
    string::size_type endt;
    // is this token a [, matlab bracket
    if ( strc.substr(idx,1) == chL ) {

     // get position of ending bracket
      if ( (idx = strc.find_first_of(chR,idx+1,lgt)) == string::npos)  {
        cerr << "VG_Pilot matlab token: no right bracket found" << endl;
        cerr << "check pilot file for matching [] brackets" << endl;
        cerr << "pilotfile line with problem follows" << endl;
        cerr << strc << endl;
        exit(0);
      }
      
      endt = idx;
      idx++;
      if ( (strc.compare(idx,1,chB,1) !=0) &&
           (idx < strc.size() ) ){
        cerr << "non-blank character following ] in matlab format" << endl;
        cerr << "check pilot file, ] should be end of token" << endl;
        cerr << "pilotfile line with problem follows" << endl;
        cerr << strc << endl;
        exit(0);
      }
    }
    else {
      // get position of next blank character
      idx = strc.find_first_of(chB,idx+1,lgt);
      endt = idx-1;
    }

    string tokenstr = strc.substr(start,endt-start+1);

    if ( (tokenstr.find("]") != string::npos) &&
         (tokenstr.find("[") == string::npos) ) {
 
      cerr << "found pilot line with ] bracket and no matchine [" << endl;
      cerr << "pilotfile line with problem follows:" << endl;
      cerr << strc << endl;
      exit(0);
    }

   tokens.push_back(tokenstr);
  }

  return; 
}
//************ end of tokenizer ******************************
