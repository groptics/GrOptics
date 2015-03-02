/*
VERSION3.1
2March2015
*/
/*! \brief GPilot class for reading pilotfiles.

GPilot provides a standard method for reading the pilot files 
used in many Grinnell/ISU programs.

*/
#ifndef GPILOT
#define GPILOT

class GPilot {
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
  /*! GPilot constructor.
    \param filename  name of pilotfile
  */
  GPilot(const string filename);   

  /*! GPilot default constructor.
  */
  GPilot();   

  /*! addPilotFile: provides filename for GPilot to read
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

