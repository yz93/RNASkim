#include <vector>
#include <map>
#include "rs_selected_key.h"
//#include "proto/rnasigs.pb.h"

using std::vector;
using std::map;
 extern int FLAGS_read_length = 40;
 extern int FLAGS_rs_length = 40;
namespace rs {
  //extern const int FLAGS_rs_length = 40;
  class SignatureInfo {
  public:
    // Total number of counts in the read data
    int total_counts;
    // The number of current k-mer class, e.g., if two k-mers both
    // show up in the same set of transcripts, this should be two.
    int occurences;
  };

  typedef map<vector<double>, SignatureInfo> SignatureInfoDB;

  void EM(int num_tids, const SignatureInfoDB& db,
          vector<double>* count_per_tid, vector<double>* pi);

  bool prepare_SignatureInfoDB(const SelectedKey &sk, SignatureInfoDB *db);
} // namesapce
