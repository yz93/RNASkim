#include "rs_selected_key.h"
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
using std::cout;
using std::endl;
//using std::fstream;
//using std::iostream;
using std::ios;
using std::map;
using std::string;
using std::vector;
using std::stoi;

namespace rs {
	std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
		std::stringstream ss(s);
		std::string item;
		while (std::getline(ss, item, delim)) {
			elems.push_back(item);
		}
		return elems;
	}

	std::vector<std::string> split(const std::string &s, char delim) {
		std::vector<std::string> elems;
		split(s, delim, elems);
		return elems;
	}


	map<string, SelectedKey> process_selected_keys(string cluster_file, string sigmer_count_file, string sigmer_details)
	{
		map<string, SelectedKey> ret;
		string line;
		string gene;
		string tid;
		string tidx;  // transcript index
		string sigmer;
		string indices;
		string sigmer_cnt;
		string length_str;
		
		std::ifstream cluster_f(cluster_file); // "Z:\\RNASkim\\RNASkim\\src\\emfiles\\cluster_info.txt"
		std::ifstream sig_cnt_f(sigmer_count_file); //"Z:\\RNASkim\\RNASkim\\src\\emfiles\\sigmer_occurance.txt"
		std::ifstream sig_det_f(sigmer_details); //"Z:\\RNASkim\\RNASkim\\src\\emfiles\\sigmer_detail.txt"

		while (getline(cluster_f, line)){
			std::istringstream iss(line);
			iss >> gene;
			iss >> tid;
			iss >> length_str;  // skip transcript index
			iss >> length_str;
			ret[gene].tids.push_back(tid);
			ret[gene].lengths.push_back(std::stoi(length_str));
			ret[gene].gid = "gene";
		}

		cluster_f.close();

		while (getline(sig_cnt_f, line)){
			std::istringstream iss(line);
			iss >> gene;
			iss >> sigmer;  // skip a random binary string
			iss >> sigmer;  
			iss >> sigmer_cnt;
			ret[gene].keys.push_back(Key(sigmer, std::stoi(sigmer_cnt)));
		}

		while (getline(sig_det_f, line)){
			std::istringstream iss(line);
			iss >> gene;
			iss >> sigmer;  // skip a random binary string and read the sigmer
			iss >> tidx;
			iss >> indices;
			indices = indices.substr(1, indices.size()-2);
			vector<string> positions_str;
			positions_str = split(indices, ',');
			// Assumes that stoi strips away the space at the beginning of the string
			vector<int> positions;
			for (size_t i = 0; i < positions_str.size(); ++i){
				//cout << stoi(positions_str[i]) << endl;
				positions.push_back(std::stoi(positions_str[i]));
			}
			// ================================================================
			TranscriptInfo temp_ti(stoi(tidx));  // create a TranscriptInfo object
			temp_ti.positions = positions;    // Assign the position vector to its attribute
			size_t i;
			for (i = 0; i < ret[gene].keys.size(); ++i){
				if (ret[gene].keys[i].key == sigmer){
					ret[gene].keys[i].transcript_infos.push_back(temp_ti);
					break;
				}
			}
			if (i == ret[gene].keys.size()){
				ret[gene].keys.push_back(Key(sigmer, 0));
			}
		}
		// By now, all SelectedKey objects should have been successfully stored in ret.
		return ret;
	}
}
/*
int main(){
	rs::SelectedKey sk;
	map<std::string, rs::SelectedKey> complete_gene_to_sk;
	string path1 = "Z:\\RNASkim\\RNASkim\\src\\emfiles\\cluster_info.txt";
	string path2 = "Z:\\RNASkim\\RNASkim\\src\\emfiles\\sigmer_detail.txt";
	string path3 = "Z:\\RNASkim\\RNASkim\\src\\emfiles\\sigmer_occurance.txt";
	complete_gene_to_sk = rs::process_selected_keys(path1, path2, path3);
	int x = 0;
	for (auto iter : complete_gene_to_sk){
		if (x == 1){
			cout << "The gene is " << iter.first << endl;
			cout << "Its corresponding transcripts are: " << endl;
			for (auto t : iter.second.tids){
				cout << t << "  ";
			}
			cout << endl;
			cout << "The transcripts corresponding lengths are: " << endl;
			for (auto t : iter.second.lengths){
				cout << t << "  ";
			}
			cout << endl;
			cout << "Its sigmers are: " << endl;
			for (auto t : iter.second.keys){
				cout << t.key << "  ";
			}
			cout << endl;
		}
		++x;
	}
}
*/
