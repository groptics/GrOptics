// version: 6.1.0    date: 2014-05-06
/*! \brief VG_Pilot class for reading pilotfiles.

VG_Pilot provides a standard method for reading the pilot files 
used in many Grinnell/ISU programs.

*/
#ifndef VG_PILOT
#define VG_PILOT

extern ostream *oLog;

class VG_Pilot {
 private:
  const string pilot_file;  //!< name of pilot file
  vector<string> tokens;    //!< container for record parameters
  vector<string> flagline;  //!< flagged record read from pilotfile

  multimap<const string,vector<string> > m; //!< map with pilot records
  multimap<const string,vector<string> >::iterator iter;
  multimap<const string,vector<string> >::iterator lower;
  multimap<const string,vector<string> >::iterator upper;
  int count;          
  int line_count;

  typedef pair<string,vector<string> > entry;
  string flag;      //!< flag following asterisk in pilotfile

  /*! tokenizer: breaks string into a vector of string tokens
    \param str: string to convert to tokens
    \param tokens: string vector to hold tokens array
  */
  void tokenizer(string& str, vector<string>& tokens);

 public:
  /*! VG_Pilot constructor.
    \param filename  name of pilotfile
  */
  VG_Pilot(const string filename);   

  /*! VG_Pilot default constructor.
  */
  VG_Pilot();   

  /*! addPilotFile: provides filename for VG_Pilot to read
    and append to its map.
   */
  bool addPilotFile(const string &pilotfile);

  /*! set_flag: sets flagname for flag search 
    \param flagname  name of flag in pilotfile
    \return number of such flag lines in pilotfile
   */
  int set_flag(const string flagname);

  /*! returns next line with flag from set_flag
    \param line_vector  vector of values following flag in pilot line
    \return remaining number of these flag lines to be read
   */
  int get_line_vector(vector<string> &line_vector);
};

#endif


/* first use the set_flag method 
   then, can use the get_line_vector method until count 
   reaches zero.  Note that count is reset following any set_flag
   call and you have to start from the beginning in counting the number
   of flag lines to read

   if flag doesn't exist, then get_line_vector returns an empty vector
   following set_flag(flag that doesn't exist).

   see test_pilot.cpp for example use with test.pilot file

*/

