#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include <deque>
#include <list>
#include <iterator>
#include <algorithm>
#include <bitset>
#include <iomanip>
#include <ctime>

using namespace std;
#include "RConfig.h"
int main(int argc, char *argv[]) {

  int *i = new int(1);

  cout << "*i = " << *i << endl;
  SafeDelete(i);
  SafeDelete(i);
}
