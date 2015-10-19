// TODO(zzj): make this program multithread

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

//#include "gflags/gflags.h"
//#include "glog/logging.h"

//#include "rs_common.h"
#include "rs_estimate_lib.h"
//#include "proto/rnasigs.pb.h"
//#include "proto_data.h"

using std::cout;
using std::endl;
using std::fstream;
using std::ios;
using std::map;
using std::string;
using std::vector;

//string path1 = "Z:\\RNASkim\\RNASkim\\src\\emfiles\\cluster_info.txt";
//string path2 = "Z:\\RNASkim\\RNASkim\\src\\emfiles\\sigmer_detail.txt";
//string path3 = "Z:\\RNASkim\\RNASkim\\src\\emfiles\\sigmer_occurance.txt";
/*
DEFINE_string(count_file, "",
              "The path to the file that contains all selected keys with their counts (output).");
DEFINE_int32(rs_length, 40,
           "The length of the sig-mer.");
DEFINE_int32(read_length, 100,
           "The length of the RNA-seq reads.");
*/

//int FLAGS_rs_length = 40;

#define DEBUG 0

namespace rs {
  class EstimateMain {
  public:
		EstimateMain(const string path1, const string path2, const string path3)
			: m_path1(path1), m_path2(path2), m_path3(path3){};

    void run() {
      SelectedKey sk;
      map<string, double> profile;
      map<string, int> tid2length;
			map<std::string, rs::SelectedKey> complete_gene_to_sk;
			complete_gene_to_sk = rs::process_selected_keys(m_path1, m_path2, m_path3);
      //while(load_protobuf_data(&istream, &sk, buffer, buffer_size)) 
			for (auto gene_iterator : complete_gene_to_sk) {
        // a table from a transcript id to a vector of estimated
        // abundunce values.
				sk = gene_iterator.second;
        bool debug = false;
        int inspected_tid = 0;
        for (int i = 0; i < sk.tids_size(); i++) {
          tid2length[sk.tids[i]] = sk.lengths[i];
          // if (sk.tids(i) == "ENSMUST00000129709") {
          if (sk.tids[i] == "ENSMUST00000114890") {
            debug = true;
            inspected_tid = i;
          }
        }
        if (DEBUG && !debug) continue;
        //LOG_IF(ERROR, sk.tids_size() > 100) << "Preparing";

        SignatureInfoDB db;
        bool run_em = prepare_SignatureInfoDB(sk, &db);

        vector<vector<int> > covered_transcripts; // = find_covered_transcript(db, sk.tids_size());
        vector<double> density_per_tid(sk.tids_size());
        //LOG_IF(ERROR, sk.tids_size() > 100) << "Start EM";
        if (run_em) {
          vector<double> pi(sk.tids_size(), 1.0 / sk.tids_size());
          EM(sk.tids_size(), db, &density_per_tid, &pi);
          // bool penalized_run = false;
          // vector<double>  penalty(sk.tids_size(), 1);
          // for (int i = 0; i < (int) covered_transcripts.size(); i++) {
          //   if (covered_transcripts[i].size() == 0) continue;
          //   if (density_per_tid[i] > 0) {
          //     penalized_run = true;
          //     for (int j : covered_transcripts[i]) {
          //       penalty[j] *= 5;
          //     }
          //   }
          // }
          // if (penalized_run) {
          //   SignatureInfoDB new_db;
          //   for (const auto& item : db) {
          //     vector<double> w = item.first;
          //     for (int i= 0; i < sk.tids_size(); i++) {
          //       w[i] *= penalty[i];
          //     }
          //     new_db[w] = item.second;
          //   }
          //   density_per_tid = EM(sk.tids_size(), new_db);
          // }
        }
        vector<double> num_reads(sk.tids_size(), 0);
        for (int j = 0; j < sk.tids_size(); j++) {
          num_reads[j] = density_per_tid[j] * (sk.lengths[j]);
        }
        for (int j = 0; j < sk.tids_size(); j++) {
          if (num_reads[j] < 0.1) num_reads[j] = 0;
        }
        if (DEBUG) {
          for (size_t i = 0; i <  covered_transcripts.size(); i++) {
            if (covered_transcripts[i].size() == 0) continue;
            cout << sk.tids[i] << " covers ";
              for (size_t j = 0; j < covered_transcripts[i].size(); j++) {
                cout << sk.tids[covered_transcripts[i][j]] << ' ';
              }
              cout << endl;
          }
        }
        for (size_t j = 0; j < num_reads.size(); j++) {
          profile[sk.tids[j]] = num_reads[j];
        }
      }
      double estimated_total_reads = 0;
      double estimated_total_abundance = 0;
      for (auto& iter : profile) {
        iter.second /= FLAGS_read_length;
        iter.second =
            iter.second * FLAGS_read_length / (FLAGS_read_length - FLAGS_rs_length);
        estimated_total_reads += iter.second;
        estimated_total_abundance += iter.second / tid2length[iter.first];
      }
      //LOG(ERROR) << "Estimated total reads " << estimated_total_reads;
      for (auto iter : profile) {
        double rpkm = 0;
        double tpm = 0;
        int length = 0;
        if (tid2length.find(iter.first) != tid2length.end()) {
          rpkm = iter.second / (tid2length[iter.first] / 1000.0) / (estimated_total_reads / 1000000.0);
          // the sum of tpm is 1 million
          tpm = iter.second / tid2length[iter.first] / estimated_total_abundance * 1000000;
          length = tid2length[iter.first];
        }
        std::cout << iter.first << '\t' << length << '\t'
                  << iter.second << '\t' << rpkm << '\t' << tpm << std::endl;
        // for (auto v : iter.second) {
        //   cout << v << '\t';
        // }
        // cout << endl;
      }
    }
  private:
    string m_path1;
		string m_path2;
		string m_path3;
  };
}  // namespace rs


int main(int argc, char *argv[]) {
  //google::InitGoogleLogging(argv[0]);
  //google::ParseCommandLineFlags(&argc, &argv, true);
	if (argc < 7) {
		std::cerr << "Usage: " << argv[0] << "--cluster_info_file PATH --sigmer_count_file PATH --sigmer_details_file PATH" << std::endl;
		return 1;
	}

	std::vector <std::string> sources;
	std::string path1;
	std::string path2;
	std::string path3;
	for (int i = 1; i < argc; ++i) {
		if (std::string(argv[i]) == "--cluster_info_file") {
			if (i + 1 < argc) { // Make sure we aren't at the end of argv!
				++i;
				path1 = argv[i];
			}
			else { 
				std::cerr << "--cluster_info_file option requires one argument." << std::endl;
				return 1;
			}
		}
		else if (std::string(argv[i]) == "--sigmer_count_file"){
			if (i + 1 < argc) { // Make sure we aren't at the end of argv!
				++i;
				path2 = argv[i];
			}
			else { // Uh-oh, there was no argument to the destination option.
				std::cerr << "--sigmer_count_file option requires one argument." << std::endl;
				return 1;
			}
		}
		else if (std::string(argv[i]) == "--sigmer_details_file"){
			if (i + 1 < argc) { // Make sure we aren't at the end of argv!
				++i;
				path3 = argv[i];
			}
			else {
				std::cerr << "--sigmer_details_file option requires one argument." << std::endl;
				return 1;
			}
		}
		else {
			std::cerr << "Wrong Format" << std::endl;
		}
	}
  rs::EstimateMain em(path1, path2, path3);
  em.run();
}
